#include "RandomGenerator.h"

RandomGenerator::RandomGenerator()
{
    std::random_device rd;
    m_rng.seed(rd());
}

RandomGenerator::~RandomGenerator() {}