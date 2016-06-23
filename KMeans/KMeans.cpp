// KMeans.cpp : Defines the entry point for the console application.
//


#include <iostream>
#include "kmeans.h"

using namespace std;
using namespace kmc;


int main(int argc,     // Number of strings in array argv
	char *argv[],      // Array of command-line argument strings
	char *envp[])      // Array of environment variable strings
{

	auto PrintHelp = [&]()
	{
		cout << "Assign 3d points to K clusters using K-Means\n\n";
		cout << "Usage: "<< argv[0] << " path_to_file K [Iter] [Tol]\n\n";
		cout << "K : number of clusters\n";
		cout << "Iter : Max number of iterations [default=100]\n";
		cout << "Tol  : Tollerance in center movement for convergence [default=0.1]\n\n";
		cout << argv[0] << " --help shows this help";
	};


	if (argc <= 2 || argc > 5)
	{
		PrintHelp();
		exit(0);
	}

	auto Data = LoadPointCloud(argv[1]);

	int K = 0;
	int MaxIter = 100;
	float Tol = 0.1 ;

	try
	{
		K = std::stoi(argv[2]);
	}
	catch (...)
	{
		std::cout << "\nERROR: Invalid K" << "\n\n";
		PrintHelp();
		exit(-1);
	}

	if (argc > 3)
	{
		try
		{
			MaxIter = std::stoi(argv[3]);
		}
		catch (...)
		{
			std::cout << "\nERROR: Invalid Iter" << "\n\n";
			PrintHelp();
			exit(-1);
		}
	}

	if (argc == 5)
	{
		try
		{
			Tol = std::stof(argv[4]);
		}
		catch (...)
		{
			std::cout << "\nERROR: Invalid Tol" << "\n\n";
			PrintHelp();
			exit(-1);
		}
	}



	if (K < 1)
	{
		std::cout << "ERROR: K must be larger than 0" << "\n";
		system("PAUSE");
		exit(-1);
	}

	if (MaxIter < 1)
	{
		std::cout << "ERROR: Iter must be larger than 0" << "\n";
		system("PAUSE");
		exit(-1);
	}

	if (K >= Data.NPoints)
	{
		std::cout << "ERROR: The number of datapoints must be larger or equal to the number of clusters" << "\n";
		system("PAUSE");
		exit(-1);
	}



  KMeans<EuclidianDistance> km(std::move(Data), K);

  if (!km)
  {
	  std::cout << "ERROR: Cannot initialize KMeans class (allocation or invalid parameters)" << "\n";
	  exit(-1);
  }

  km.Cluster(MaxIter,Tol);
  km.PrintCenters();
  system("PAUSE");
  
}

