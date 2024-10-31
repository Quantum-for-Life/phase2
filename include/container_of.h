#ifndef CONTAINER_OF_H
#define CONTAINER_OF_H

#define container_of(ptr, type, membr)                                         \
	((type *)((void *)(ptr) - offsetof(type, membr)))

#endif // CONTAINER_OF_H
