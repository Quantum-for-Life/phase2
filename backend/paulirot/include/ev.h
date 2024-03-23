#ifndef EV_H
#define EV_H

struct ev {
	int num_ranks;
	int rank;
};

int ev_init(struct ev *);

int ev_destroy(struct ev *);

#endif // EV_H
