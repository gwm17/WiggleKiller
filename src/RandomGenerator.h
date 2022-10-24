#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <random>

class RandomGenerator
{
public:
    ~RandomGenerator();

    std::mt19937& GetGenerator() { return m_rng; }
    static RandomGenerator& GetInstance()
    {
        static thread_local RandomGenerator generator;
        return generator;
    }

private:
    RandomGenerator();

    std::mt19937 m_rng;
};

#endif