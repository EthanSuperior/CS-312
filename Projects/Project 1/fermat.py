import random
from math import floor

def prime_test(N, k):
	# This is main function, that is connected to the Test button. You don't need to touch it.
	return fermat(N,k), miller_rabin(N,k)

# This recursive function finds the modular exponential of a number in the form: x^y % N
# We do this by the realization that exponential can be written as the sum of their parts
# ie x^12 is the same as x^10 + x^2 or, as it is more useful x^8 + x^4
# by taking a given value and doubling it multiple times we can quickly get x^y
# by splitting the number where appropriate we can than add in the sum of the parts by multiplying by x
# anytime we would have added, in the case of the x^12 we can write it as x^8 + x^4
# this algorithm would evaluate this as (x * (1 ^2))) ^2 ^2 ^2 which is x^8, and to deal with the missing x^4
# We would multiply by an additional x in the midst of the doubling
# Yielding (x * ((x * (1 ^2)) ^2)) ^2 ^2 conveniently if we half the value to which we are trying to raise
# x to and then check if it is odd or even, we get the necessary times of which to multiply in this additional x

# We use z to represent the already evaluated recursive segments
# Written recursively it evaluates backwards as follows: x^12 = (x * ((x * (1 ^2)) ^2)) ^2 ^2
#  CALL            RETURN
# x, y = 12     z (x^6) is squared and as y is even we return it as is, ending the recursion at x^12
# x, y = 6      z (x^3) is squared, y is even, so it is returned as is x^6
# x, y = 3      z (x) is squared, again the value of y is odd, so we add in another x, this is x^3 (x*((x*(1^2))^2))
# x, y = 1      The inner z (1), is returned and squared, y is odd, so we multiply by x, getting x as the returning z
# x, y = 0      Anything to the zeroth power is one, hence the return 1 closing the recursion and returning z = 1

# By taking the modulus of these numbers, we can find the modular exponential

# Time: O(m^3) Space: O(m^2)
def mod_exp(x, y, N):
    # m calls, each taking 3m space -> O(m^2)
    # For the base case where x^0 is always 1
    if y == 0:  # O(1)
        return 1
    # Recursive call with a reduction to y, driving it towards the base case
    z = mod_exp(x, floor(y / 2), N)  # mod_exp: O(m^#), Div: O(m^2), Calls: m times
    # If we use bitwise AND we can quickly check if the number is odd or even by looking at the final bit
    # the ones bit, if it has a 1, then the bitwise AND will return 1, which evaluates to true
    # So if the number is odd (1 in the final binary bit) it will multiply
    # in the missing x needed to account for the addition of the next successive power [the x^4 in the above example]
    if y & 1:
        return (x * (z * z)) % N  # Mod: O(m^2), Mult: O(m^2) Space: O(m), Mult O(m)
    # Otherwise, as the number is even, we can know with surety that the z contains all the needed x's
    else:
        return (z * z) % N  # Mod: O(m^2), Mult: O(m^2)
	

# The probability that the Fermat Primality works on a given a is 50%
# So for k tests of a you get a further 50% reduction
# Therefore 1 - .5^k will give you your probability of it actually being a prime
# Time: O(k^3) Space: O(k^2)
def fprobability(k):
    return 1.0 - mod_exp(.5, k, k)  # mod_exp: O(m^3), Sub: O(m); Space: mod_exp: O(m^2), Sub O(m)

# The probability that the Miller-Rabin Primality works on a given a is 75%
# So for k tests of a you get a further 75% reduction
# Therefore 1 - .75^k will give you your probability of it actually being a prime
# Time: O(k^3) Space: O(k^2)
def mprobability(k):
    return 1.0 - mod_exp(.25, k, k)  # mod_exp: O(m^3), Space: O(m^2)

# Fermat's Theorem states that for a given number N
# we can know that a number is NOT prime when a value a, such that 1<=a<N,
# a^(N-1) % N does not equal 1
# conversely it has a high probability (at least 50%) of being prime if it does

# We implemented this theorem with k trials to increase the probability of success
# Time: O(m^3*log2(m)) Space: O(m^2)
def fermat(N,k):
    for i in range(k):  # O(k*(m^3)) k is at most log2(n), where m is the length of N in binary
        a = random.randint(1, N - 1)  # O(1)
        if mod_exp(a, N - 1, N) != 1:  # mod_exp: O(m^3), space is O(m^2)
            return 'composite'  # O(1)
    return 'prime'  # O(1)

# Miller-Rabin's Test improves upon the Fermat's Theorem
# It does this by further testing a passing a value,
# If you take the square root of that a, it should also pass the test,
# if the a does pass, you continue taking the square root and testing it until eventually it fails.
# Upon failing you check to see if the new modulus remaining is only 1 short, ie N-1
# if it is than the number is prime with a 75% accuracy, otherwise it is proven to actually be composite

# We implemented this test with k trials to increase the probability of success
# Time: O(m^4log2(m)) Space: O(m^2)
def miller_rabin(N,k):

    for i in range(k):  # O(k*(m^3)) k is at most log2(n), where m is the length of N in binary
        a = random.randint(1, N - 1)  # O(1)
        mod = mod_exp(a, N - 1, N)  # O(m^3)
        if mod != 1:  # O(1)
            return 'composite'  # O(1)
        else:
            b = N - 1  # O(m)
            while mod == 1 and not b & 1:  # Max m times
                b = b // 2   # O(m^2)
                mod = mod_exp(a, b, N)  # O(m^3), Space:O(m^2)
            if mod == N - 1 or mod == 1:  # O(1)
                return 'prime'  # O(1)
    return 'composite'  # O(1)
