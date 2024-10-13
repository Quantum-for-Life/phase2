#ifndef WORLD_CUDA_H
#define WORLD_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif


struct world_cuQuantum {
	int local_rank;
	int local_size;
};


#ifdef __cplusplus
}
#endif

#endif /* WORLD_CUDA_H */
