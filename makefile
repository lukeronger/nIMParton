##### IMParton v0.0
##### Authors : Rong Wang, Baiyang Zhang, Yun Xue
##### date : 2013-09-23
##### IMParton v1.0
##### Authors : Rong Wang
##### updated at 2016-07-17
##### contact : rwang@impcas.ac.cn / rwangcn8@gmail.com
##### nIMParton v1.0
##### updated at 2016-10-07
##### Author : Rong Wang,  rwang@impcas.ac.cn / rwangcn8@gmail.com

all : test
test : nIMParton.o test.o
	g++ -Wall -s -o $@ nIMParton.o test.o
%.o : %.cpp %.h
	g++ -O2 -fPIC -DLINUX -Wall -c -o $@ $<

.PHONY : clean
clean:
	rm -f test nIMParton.o test.o
