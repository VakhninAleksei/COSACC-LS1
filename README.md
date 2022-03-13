# Coordination of Self-adaptive Cooperative Coevolution algorithms with Local Search
This program code implements the Coordination of Self-adaptive Cooperative Coevolution algorithms with Local Search (COSACC-LS1) algorithm for solving large-scale global optimization problems. The proposed self-adaptation method tunes both the structure of the complete approach and the parameters of each algorithm in the cooperation. The SHADE algorithm is used as a subcomponent optimizer of the CC-based algorithms. Multiple Trajectory Search (MTS-LS1) uses as a local search algorithm. The problems are presented as the benchmark set (LSGO CEC'2013). We use the MPICH2 framework to implement parallel numerical experiments on our computing cluster.

# MPICH2
This programming code works using the MPICH2 framework. Information about intallation of MPICH2 can be found in this [web site](https://mpitutorial.com/tutorials/installing-mpich2/). The "hostfile" contains information about PC-slaves in out local network, thair IP-addresses and how many computational thread we will use. 


```
192.168.1.110:16 #master
192.168.1.35:16 #slave1
192.168.1.39:16 #slave2
192.168.1.45:16 #slave3
192.168.1.99:16 #slave4
192.168.1.50:16 #slave5
192.168.1.81:16 #slave6
192.168.1.118:16 #slave7
```
In our case, we have eight PC

We use "perform.sh" file to run compiled program.
```
mpirun -np 128 --hostfile hostfile ./COSACC-LS1
```
