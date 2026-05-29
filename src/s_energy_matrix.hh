#ifndef ENERGY_MATRIX_H
#define ENERGY_MATRIX_H

#include "base_types.hh"
#include "matrices.hh"
#include "hotspot.hh"
#include "SHAPE.hh"
#include <string>
#include <vector>

extern "C" {
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params/io.h"
}
bool compare_hotspot_ptr(Hotspot &a, Hotspot &b);
bool compare_hotspot_ptr_pos(Hotspot &a, Hotspot &b);

class s_energy_matrix {
  public:

    s_energy_matrix(std::string &seq, SHAPEData *ShapeData);
    // The constructor

    ~s_energy_matrix();
    // The destructor

    vrna_param_t *params_;

    short *S_;
    short *S1_;
    SHAPEData *ShapeData;

    // return the type at V(i,j)
    // Mateo 13 Sept 2023
    void compute_hotspot_energy(cand_pos_t i, cand_pos_t j, bool is_stack);

    energy_t compute_stack(cand_pos_t i, cand_pos_t j, const paramT *params);
    energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params);
    void get_hotspots(std::vector<Hotspot> &hotspot_list, int max_hotspot, int min_stem, int turn);
    void expand_hotspot(Hotspot &hotspot, int n);

    // better to have protected variable rather than private, it's necessary for Hfold
  protected:
    // private:

    std::string seq_;
    cand_pos_t n; // sequence length
    std::vector<cand_pos_t> index;
    TriangleMatrix energy;
};

#endif
