#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H


struct rayon_circuit_data {
        PauliHamil hamil;
        void *data; // state preparation. TBA
};

struct rayon_circ_data {
        double time;
        int imag_switch;
};

struct circuit rayon_circuit_factory(struct rayon_circuit_data *);

#endif //PHASE2_RAYON_H
