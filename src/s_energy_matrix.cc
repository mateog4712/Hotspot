/***************************************************************************
                          s_energy_matrix.cpp  -  description
                             -------------------
    begin                : Fri Apr 12 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#include "s_energy_matrix.hh"

s_energy_matrix::s_energy_matrix(std::string &seq,  SHAPEData *ShapeData) : params_(scale_parameters())
// The constructor
{
    make_pair_matrix();
    S_ = encode_sequence(seq.c_str(), 0);
    S1_ = encode_sequence(seq.c_str(), 1);

    n = seq.length();
    seq_ = seq;
    this->ShapeData = ShapeData;

    energy.new_index(index,n+1);
    energy.init(n+1,index);
}

s_energy_matrix::~s_energy_matrix()
// The destructor
{
    free(params_);
    free(S_);
    free(S1_);
}

energy_t s_energy_matrix::compute_stack(cand_pos_t i, cand_pos_t j, const paramT *params) {

    const int ptype_closing = pair[S_[i]][S_[j]];
    cand_pos_t k = i + 1;
    cand_pos_t l = j - 1;
    return E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1],
                     const_cast<paramT *>(params))
           + energy.get(k, l) + ShapeData->get_calculated(i) + ShapeData->get_calculated(j);
}

energy_t s_energy_matrix::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params) {

    const int ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1],
                     const_cast<paramT *>(params))
           + energy.get(k, l);
}

// Mateo May 2026
void s_energy_matrix::compute_hotspot_energy(cand_pos_t i, cand_pos_t j, bool is_stack) {
    energy_t en = 0;
    if (is_stack) {
        en = compute_stack(i, j, params_);
        en+=ShapeData->get_calculated(i);
        en+=ShapeData->get_calculated(j);
    } else {
        en = 0; // HairpinE(seq_,S_,S1_,params_,i,j);
    }
    energy.set(i,j) = en;
    return;
}

// Mateo May 2026
// given a initial hotspot which is a hairpin loop, keep trying to add a arc to form a larger stack
void s_energy_matrix::expand_hotspot(Hotspot &hotspot, int n) {
    // calculation for the hairpin that is already in there
    compute_hotspot_energy(hotspot.get_left_outer_index(), hotspot.get_right_outer_index(), 0);

    // try to expand by adding a arc right beside the current out most arc
    while (hotspot.get_left_outer_index() - 1 >= 1 && hotspot.get_right_outer_index() + 1 <= n) {
        base_type sim1 = S_[hotspot.get_left_outer_index() - 1];
        base_type sjp1 = S_[hotspot.get_right_outer_index() + 1];
        pair_type ptype_closing = pair[sim1][sjp1];
        if (ptype_closing > 0) {
            hotspot.move_left_outer_index();
            hotspot.move_right_outer_index();
            hotspot.increment_size();
            compute_hotspot_energy(hotspot.get_left_outer_index(), hotspot.get_right_outer_index(), 1);
        } else {
            break;
        }
    }
    base_type i = hotspot.get_left_outer_index();
    base_type j = hotspot.get_right_outer_index();
    pair_type tt = pair[S_[i]][S_[j]];
    base_type si1 = i > 1 ? S_[i - 1] : -1;
    base_type sj1 = j < n ? S_[j + 1] : -1;
    energy_t dangle_penalty = vrna_E_ext_stem(tt, si1, sj1, params_);

    double en = energy.get(hotspot.get_left_outer_index(), hotspot.get_right_outer_index());

    en = (en + dangle_penalty) / 100;

    hotspot.set_energy(en);
    return;
}

// Mateo May 2026
// look for every possible hairpin loop, and try to add a arc to form a larger stack with at least min_stack_size bases
void s_energy_matrix::get_hotspots(std::vector<Hotspot> &hotspot_list, int max_hotspot, int min_stem, int turn) {

    int min_stack_size = min_stem; // the hotspot must be a stack of size >= min_Stem
    // start at min_stack_size-1 and go outward to try to add more arcs to form bigger stack because we cannot expand more than min_stack_size from
    // there anyway
    int max_i = n-(min_stack_size+turn)+1;
    int min_j = n-min_stack_size+1;
    for (cand_pos_t i = min_stack_size; i <= max_i; i++) {
        for (cand_pos_t j = i+turn+1; j <= min_j; j++) {
            int ptype_closing = pair[S_[i]][S_[j]];
            if (ptype_closing > 0) {

                Hotspot current_hotspot(i, j, n);

                expand_hotspot(current_hotspot, n);

                if (current_hotspot.get_size() < min_stack_size || current_hotspot.is_invalid_energy()) {

                } else {

                    current_hotspot.set_structure();
                    hotspot_list.push_back(current_hotspot);
                }
            }
        }
    }

    // make sure we only keep top max_hotspot hotspots with lowest energy
    std::sort(hotspot_list.begin(), hotspot_list.end(), compare_hotspot_ptr);
    while ((int)hotspot_list.size() > max_hotspot) {
        hotspot_list.pop_back();
    }

    return;
}

bool compare_hotspot_ptr(Hotspot &a, Hotspot &b) { return (a.get_energy() < b.get_energy()); }
