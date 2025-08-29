import sys
import os
import glob
import random
import math
import numpy as np
# prime number, adapted from Chatgpt
def is_prime(n):
    """Checks if a number is prime using trial division."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def random_large_prime(bit_length):
    """Generates a random prime number of the given bit length."""
    while True:
        num = random.getrandbits(bit_length) | 1  # Ensure it's odd
        if is_prime(num):
            return num

def generate_prime_list(count, bit_length, minValue):
    """Generates a list of 'count' random large prime numbers."""
    primes = []
    while len(primes) < count:
        prime = random_large_prime(bit_length)
        if prime > minValue:
            if prime not in primes:  # Ensure uniqueness
                primes.append(prime)
    return primes

def computeSphereVolume(sphereRadius):
    return (4/3) * np.pi * (sphereRadius **3 )