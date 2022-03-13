# Coordination of Self-adaptive Cooperative Coevolution algorithms with Local Search (COSACC-LS1)
This program code implements the Coordination of Self-adaptive Cooperative Coevolution algorithms with Local Search (COSACC-LS1) algorithm for solving large-scale global optimization problems. The proposed self-adaptation method tunes both the structure of the complete approach and the parameters of each algorithm in the cooperation. The SHADE algorithm is used as a subcomponent optimizer of the CC-based algorithms. Multiple Trajectory Search (MTS-LS1) uses as a local search algorithm. The problems are presented as the benchmark set (LSGO CEC'2013). We use the MPICH2 framework to implement parallel numerical experiments on our computing cluster.

# COSACC-LS1 settings
The main settings of COSACC-LS1 are:
```
generations_init = 20; // the number of initial generations for all CC-based algorithms
pop_size_min = 25; // the minimum number of individuals
pop_size_init = 100; // the initial number of individuals
pop_size_max = 200; // the maximum number of individuals
FEV_LS1_budget = 25000; // the number of fitness evaluations for MTS-LS1
island_setup [8][3] = // different settings of COSACC for numerical experiment
    {1,2,4}, {1,2,8}, {1,2,10}, {1,4,8}, {1,4,10}, {1,8,10}, {2,4,8}, {2,4,10} // the sets of CC-based algorithms
};
```

# MPICH2
This programming code works using the MPICH2 framework. Information about installation of MPICH2 can be found in this [web site](https://mpitutorial.com/tutorials/installing-mpich2/). The "hostfile" contains information about PC-slaves in our local network, their IP-addresses and the number of computational threads.

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

We use "perform.sh" file to run the compiled program.
```
mpirun -np 128 --hostfile hostfile ./COSACC-LS1
```
