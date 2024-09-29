#ifndef WORLD_H
#define WORLD_H

struct world {
	int num_ranks;
	int rank;
};

int world_init(struct world *wd);

#endif /* WORLD_H */
