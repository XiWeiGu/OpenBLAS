for i in {128..256}
do
	cd /root/OpenBLAS && git checkout -- param.h
	sed -i "s/#define DGEMM_DEFAULT_R 256/#define DGEMM_DEFAULT_R $i/g" param.h
	make clean && make NO_LAPACK=1 USE_SIMPLE_THREADED_LEVEL3=1 NO_AFFINITY=0 && make install PREFIX=/opt/OpenBLAS 
	cd /root/hpl-2.2 && make arch=Linux_Loongson3a clean && make arch=Linux_Loongson3a
        cd bin/Linux_Loongson3a && ./xhpl | tee log.$i
done
