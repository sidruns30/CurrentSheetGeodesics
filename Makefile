all: metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp main.cpp
	g++-12 -o exec -O2 -std=c++17 metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp main.cpp -fopenmp

debug: metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp main.cpp
	g++-12 -o exec -O2 -std=c++17 metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp main.cpp

cnpy:metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp cnpy/cnpy.cpp main.cpp
	g++-12  -fno-common -fopenmp -o exec -O2 -std=c++17 metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp cnpy/cnpy.cpp main.cpp -fopenmp -L/usr/local/lib -lz

cnpy_deb: metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp cnpy/cnpy.cpp main.cpp
	g++-12  -fno-common -o exec -Og -std=c++17 metric/BLmetric.cpp metric/SCmetric.cpp metric/KSmetric.cpp metric/MKSmetric.cpp metric/metric.cpp fluid/getk.cpp defs.cpp geodesics/geodesic.cpp input/load_txt.cpp fluid/BHAC_MHD.cpp partition.cpp cnpy/cnpy.cpp main.cpp -L/usr/local/lib -lz
