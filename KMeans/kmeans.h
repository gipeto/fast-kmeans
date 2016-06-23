#pragma once
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>


namespace utils
{

	template<class T>
	struct aligned_deleter
	{
		void operator()(T * p) const
		{
			_aligned_free(p);
			p = nullptr;
		}
	};

	template<typename T>
	using unique_ptr_aligned = typename std::unique_ptr<T[], aligned_deleter<T>>;


	template<typename T, size_t Alignment = alignof(T)>
	auto make_unique_aligned(size_t n) -> unique_ptr_aligned<T>
	{

		auto p = _aligned_malloc(n*sizeof(T), Alignment);
		if (nullptr == p)
		{
			return unique_ptr_aligned<T>();
		}

		std::memset(p, 0, n*sizeof(T));
		return unique_ptr_aligned<T>(static_cast<T*>(p));
	}

}



namespace kmc
{

	using namespace utils;
	using namespace std;

	template<typename T>
	struct Point4
	{
		T x;
		T y;
		T z;
		T w;

		Point4 & operator +=(const Point4 & Other)
		{
			x += Other.x;
			y += Other.y;
			z += Other.z;
			return *this;
		}

		Point4 & operator -=(const Point4 & Other)
		{
			x -= Other.x;
			y -= Other.y;
			z -= Other.z;
			return *this;
		}

		Point4 & operator /=(T Scalar)
		{
			x /= Scalar;
			y /= Scalar;
			z /= Scalar;
			return *this;
		}


	};

	using Point4f = Point4<float>;
	using PointBuffer = unique_ptr_aligned<Point4f>;

	struct EuclidianDistance
	{
		inline static constexpr float Get(const Point4f & left, const Point4f & right)
		{
			return  (left.x - right.x) *  (left.x - right.x) + (left.y - right.y) *  (left.y - right.y) + (left.z - right.z) *  (left.z - right.z);
		}

	};


	struct DataPoints
	{
		unique_ptr_aligned<Point4f> Buffer;
		size_t NPoints;
	};


	auto LoadPointCloud(const std::string filename)->DataPoints
	{

		std::ifstream file(filename.c_str());

		if (!file.good())
		{
			std::cout << "ERROR: loading (" << filename << ") file not found or not good" << "\n";
			system("PAUSE");
			exit(0);
		}

		auto NLines = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
		auto PointBuffer = make_unique_aligned<Point4f, 32>(NLines);

		if (!PointBuffer)
		{
			std::cout << "ERROR: cannot allocate point cloud in memory" << "\n";
			system("PAUSE");
			exit(0);
		}

		file.clear();
		file.seekg(0, std::ios::beg);

		char buffer[255];
		auto lcnt = 0;

		while (!file.getline(buffer, 255).eof())
		{

			if (3 != sscanf_s(buffer, "%f %f %f", &PointBuffer[lcnt].x, &PointBuffer[lcnt].y, &PointBuffer[lcnt].z))
			{
				std::cout << "ERROR: unexpected file format" << "\n";
				exit(-1);
			}

			lcnt++;
		}

		return{ std::move(PointBuffer),static_cast<size_t>(NLines) };

	}


	/* 
	Implement K-Means based on :
	G. Hamerly, Making k-means even faster
	In Proceedings of the 2010 SIAM international conference on data mining (SDM 2010) (2010)  Key: citeulike:8266815 
	*/


	template<class Distance = EuclidianDistance>
	class KMeans
	{

		using T = float;
		DataPoints m_Data;
		size_t m_K = 0;

		PointBuffer m_Centers;				    // Cluster centers
		PointBuffer m_ClusterSum;			    // Sum of points in each cluster
		unique_ptr< int[]> m_ClusterPoints;     // Number of points per cluster
		unique_ptr<T[]> m_CenterMotion;		    // Keep track of the center motion between iterations
		unique_ptr<T[]> m_CenterDistance;	    // distance to the other closest center

		unique_ptr<int[]> m_Assignment;         // point membership to a cluster
		unique_ptr<T[]> m_UBound;               // upper bound on the distance between a point and its cluster center
		unique_ptr<T[]> m_LBound;               // lower bound on the distance between a point and its second closest cluster center

		std::random_device rd;
		std::mt19937 gen;
		

	public:

		KMeans(DataPoints && pts, size_t K = 3)
			: m_Data(std::move(pts)),
			m_K(K),
			m_Centers(make_unique_aligned<Point4f, 32>(K)),
			m_ClusterSum(make_unique_aligned<Point4f, 32>(K)),
			m_ClusterPoints(make_unique<int[]>(K)),
			m_CenterMotion(make_unique<T[]>(K)),
			m_CenterDistance(make_unique<T[]>(K)),
			m_Assignment(make_unique<int[]>(m_Data.NPoints)),
			m_UBound(make_unique<T[]>(m_Data.NPoints)),
			m_LBound(make_unique<T[]>(m_Data.NPoints)),
			gen(rd())
		{}


		explicit operator bool()
		{
			return (m_K > 0) && (m_K <= m_Data.NPoints)  && m_Centers && m_ClusterSum && m_ClusterPoints && m_CenterMotion && m_CenterDistance && m_Assignment && m_UBound && m_LBound;
		}



		inline void Cluster(size_t MaxIter = 20, T eps = 0.1)
		{
			
			eps = abs(eps);

			if (1 == m_K)
			{ 
				SetToAverage();
				return;
			}

			Init();
			
			auto Iter = 0;

			while (Iter < MaxIter)
			{
				UpdateCenterDistance();

				for (auto ip = 0; ip < m_Data.NPoints; ip++)
				{	
					// u < dmin1 , l > dmin2 -> u < l continue

					auto m = std::max(m_CenterDistance[ip] / 2, m_LBound[ip]);

					if (m_UBound[ip] > m)
					{
						auto pt = m_Data.Buffer[ip];

						m_UBound[ip] = Distance::Get(pt, m_Centers[m_Assignment[ip]]);

						if (m_UBound[ip] > m)
						{
							auto a_old = m_Assignment[ip];
							PointToAllCenters(pt, ip);

							if (a_old != m_Assignment[ip])
							{
								m_ClusterSum[m_Assignment[ip]] += pt;
								m_ClusterSum[a_old] -= pt;
								m_ClusterPoints[m_Assignment[ip]]++;
								m_ClusterPoints[a_old]--;	
							}
						}
					}
				}


				MoveCenters();
				UpdateBounds();

#ifdef _DEBUG
				cout << "Iteration "<< Iter <<" :\n-------------\n";
				PrintCenters();
#endif

				if (Converged(eps))
				{
					break;
				}

				Iter++;

			}

		}

		inline void PrintCenters()
		{
			cout << "\n Cluster Centers [x y z] : Number of Elements :\n";
			for (auto ic = 0; ic < m_K; ic++)
			{
				cout << "[" << m_Centers[ic].x << " " << m_Centers[ic].y << " " << m_Centers[ic].z << "] : "<< m_ClusterPoints[ic]<<"\n";
			}

		}

		private:

			inline void UpdateCenterDistance()
			{

				for (auto ic = 0; ic < m_K; ic++)
				{
					auto dmin = std::numeric_limits<T>::max();

					for (auto jc = 0; jc < m_K; jc++)
					{
						auto d = Distance::Get(m_Centers[ic], m_Centers[jc]);
						if (d < dmin && ic != jc)
						{
							m_CenterDistance[ic] = d;
							dmin = d;
						}

					}
				}


			}

			inline void PointToAllCenters(const Point4f & pt, size_t idx)
			{

				auto dmin  = Distance::Get(pt, m_Centers[0]);
				auto dmin2 = Distance::Get(pt, m_Centers[1]); 

				auto ifirstcc = 0;

				if (dmin2 < dmin)
				{
					std::swap(dmin, dmin2);
					ifirstcc = 1;
				}

				for (auto ic = 2; ic < m_K; ic++)
				{
					auto d = Distance::Get(pt, m_Centers[ic]);

					if (d < dmin)
					{
						dmin2 = dmin;
						dmin = d;
						ifirstcc = ic;
						continue;
					}

					if (d < dmin2)
					{
						dmin2 = d;
					}

				}


				m_Assignment[idx] = ifirstcc;
				m_UBound[idx] = dmin;
				m_LBound[idx] = dmin2;

			}


		

		inline void MoveCenters()
		{

			for (auto ic = 0; ic < m_K; ic++)
			{
				auto c_old = m_Centers[ic];
				auto csum  = m_ClusterSum[ic];

			// if a cluster gets empty pick a random point for the most populated cluster as new center
			if ( 0== m_ClusterPoints[ic])
			{

#ifdef _DEBUG
				cout << "Cluster " << ic << " is empty!\n";
#endif
			 	
				auto most_populated = 0;
				auto maxpts = 0;
				for (auto iic = 0; iic < m_K; iic++)
				{
					if (m_ClusterPoints[iic] > maxpts)
					{
						maxpts = m_ClusterPoints[iic];
						most_populated = iic;
					}
		
				}
				
				std::uniform_int_distribution<> dist(0, maxpts-1);
				auto nskip = dist(gen);
		
				for (auto ia = 0; ia < m_Data.NPoints; ia++)
				{
					if (m_Assignment[ia] == most_populated )
					{
						if (nskip == 0)
						{

							auto pt = m_Data.Buffer[ia];

							m_ClusterSum[most_populated] -= pt;
							m_ClusterPoints[most_populated]--;

							m_Centers[ic]    =pt;
							m_ClusterSum[ic] =pt;
							csum = pt;
							m_ClusterPoints[ic] = 1;
							m_Assignment[ia] = ic;

							break;
						
						}

						nskip--;
					}
								
				}
	
			}

				m_Centers[ic] = (csum /= static_cast<T>(m_ClusterPoints[ic]));
				m_CenterMotion[ic] = Distance::Get(m_Centers[ic], c_old);
			}


		}


		inline void UpdateBounds()
		{

			// Get the two centers that moved the most

			auto dmax  = m_CenterMotion[0];
			auto dmax2 = m_CenterMotion[1];

			auto ifirstfurthest  = 0;
			auto iseconfurthest  = 1;


			if (dmax2 > dmax)
			{
				std::swap(dmax, dmax2);
				ifirstfurthest = 1;
				iseconfurthest = 0;
			}

			for (auto ic = 2; ic < m_K; ic++)
			{
				auto d = m_CenterMotion[ic];

				if (d > dmax)
				{
					dmax2 = dmax;
					dmax = d;
					iseconfurthest = ifirstfurthest;
					ifirstfurthest = ic;
					continue;
				}

				if (d > dmax2)
				{
					dmax2 = d;
					iseconfurthest = ic;
				}

			}


			for (auto iu = 0; iu < m_Data.NPoints; iu++)
			{
			
				m_UBound[iu] += m_CenterMotion[m_Assignment[iu]];

				if (m_Assignment[iu] == ifirstfurthest)
				{
					m_LBound[iu] -= dmax2;
					continue;
				}

				m_LBound[iu] -= dmax;
			
			}
		
		}

	

		inline void InitCenters()
		{
		
			auto pps = m_Data.NPoints / m_K;

			for (auto ic = 0; ic < m_K; ic++)
			{
				std::uniform_int_distribution<> dis(ic*pps, (ic + 1)*pps - 1);
				m_Centers[ic] = m_Data.Buffer[dis(gen)];
			}
		}


		inline void Init()
		{
			InitCenters();
			
			memset(m_ClusterPoints.get(), 0, m_K * sizeof(int));
			memset(m_ClusterSum.get(), 0, m_K * sizeof(Point4f));

			for (auto ip = 0; ip < m_Data.NPoints; ip++)
			{
				auto pt = m_Data.Buffer[ip];
				PointToAllCenters(pt, ip);
				m_ClusterSum[m_Assignment[ip]] += pt;
				m_ClusterPoints[m_Assignment[ip]]++;
			}

#ifdef _DEBUG
			cout << "Initialization:\n-------------\n";
			PrintCenters();
			cout << "\n-------------\n";
#endif

		}

		inline bool Converged(T eps)
		{
		
			for (auto ic = 0; ic < m_K; ic++)
			{
				if (m_CenterMotion[ic] > eps)
				{
					return false;
				}
			}

#ifdef _DEBUG
			cout << "\nConverged!\n";
#endif
			return true;
	
		}


		inline void SetToAverage()
		{
		
			for (auto ip = 0; ip < m_Data.NPoints; ip++)
			{
				m_Centers[0] += m_Data.Buffer[ip];
			}

			m_Centers[0] /= static_cast<T>(m_Data.NPoints);
			m_ClusterPoints[0] = m_Data.NPoints;
		
		}


	

	};



}