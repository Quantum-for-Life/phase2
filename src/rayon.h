#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

typedef struct rayon_circuit_data_ {
    PauliHamil hamil;
    void *data; // state preparation. TBA
} rayon_circuit_data;

typedef struct rayon_circ_data_ {
    double time;
    int imag_switch;
} rayon_circ_data;

circuit
rayon_circuit_factory(rayon_circuit_data *);

#endif //PHASE2_RAYON_H
