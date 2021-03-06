//
// Created by elturpin on 04/12/2020.
//

#pragma once

#include <cstdint>
#include <array>

constexpr int8_t NB_BASE = 2;
constexpr int8_t CODON_SIZE = 3;
// promoter
constexpr int8_t PROM_MAX_DIFF = 4;
constexpr int8_t PROM_SIZE = 22;

const uint32_t PROM_SEQ = 0b0110100100111001101010;
const uint32_t PROM_MASK = 0b1111111111111111111111;

// terminator
constexpr int8_t TERM_STEM_SIZE = 4;
constexpr int8_t TERM_LOOP_SIZE = 3;
constexpr int8_t TERM_SIZE = TERM_STEM_SIZE + TERM_LOOP_SIZE + TERM_STEM_SIZE;

constexpr uint16_t TERM_MASK = 0b11110001111;

const std::array<uint8_t, 2048> TERM_DIST_LOOKUP = []()
{
    std::array<uint8_t, 2048> dists;

    for (size_t i = 0; i < 256; i++)
    {
        uint16_t val = i << TERM_LOOP_SIZE & 0b11110000000 | i & 0b00000001111;
        uint8_t dist = 0;

        for (size_t j = 0; j < TERM_STEM_SIZE; j++)
        {
            int right = j;
            int left = (TERM_SIZE - 1) - j;

            uint8_t left_bit = val >> right & 1;
            uint8_t right_bit = val >> left & 1;
            dist += (left_bit != right_bit);
        }
        dists[val] = dist;
    }
    return dists;
}();

// shine dalgardo
constexpr int8_t SHINE_DAL_SIZE = 6;
constexpr int8_t SD_START_SPACER = 4;
constexpr int8_t SD_TO_START = SHINE_DAL_SIZE + SD_START_SPACER + CODON_SIZE;

constexpr const uint16_t SHINE_DAL_SEQ = 0b0000000110110;
constexpr const uint16_t SHINE_DAL_SEQ_MASK = 0b1110000111111;
// stop
constexpr const uint8_t PROTEIN_END = 0b100; // CODON_STOP
constexpr const uint8_t PROTEIN_END_MASK = 0b111;

// Associate a codon sequence to its value (corresponding to its reverse binary form, eg. 110 -> 011)
const std::array<uint8_t, 8> CODON_VALUE_LOOKUP = {0b000, 0b100, 0b010, 0b110, 0b001, 0b101, 0b011, 0b111};
constexpr const uint8_t CODON_MASK = 0b111;

constexpr int32_t DO_TRANSLATION_LOOP = SHINE_DAL_SIZE + SD_START_SPACER + 3 * CODON_SIZE;

// Codon
constexpr int8_t CODON_START = 0b000;
constexpr int8_t CODON_STOP  = 0b001;
constexpr int8_t CODON_M0    = 0b100;
constexpr int8_t CODON_M1    = 0b101;
constexpr int8_t CODON_W0    = 0b010;
constexpr int8_t CODON_W1    = 0b011;
constexpr int8_t CODON_H0    = 0b110;
constexpr int8_t CODON_H1    = 0b111;

// Protein / Fuzzy space
constexpr double X_MIN = 0.0;
constexpr double X_MAX = 1.0;
constexpr double Y_MIN = 0.0;
constexpr double Y_MAX = 1.0;
constexpr double H_MIN = -1.0;
constexpr double H_MAX = 1.0;
constexpr double W_MIN = 0.0;
constexpr double W_MAX = 0.1;

constexpr int FUZZY_SAMPLING = 300;
constexpr int SELECTION_PRESSURE = 1000;

// Selection
constexpr int8_t NEIGHBORHOOD_WIDTH  = 3;
constexpr int8_t NEIGHBORHOOD_HEIGHT = 3;
constexpr int8_t NEIGHBORHOOD_SIZE   = NEIGHBORHOOD_HEIGHT * NEIGHBORHOOD_WIDTH;
