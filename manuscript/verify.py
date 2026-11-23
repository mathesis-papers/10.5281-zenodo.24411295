import argparse
import hashlib
import json
from math import comb
from pathlib import Path
from fractions import Fraction


def trim(polynomial):
    values = list(polynomial)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values or [0]


def add(left, right):
    values = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        values[index] += value
    for index, value in enumerate(right):
        values[index] += value
    return trim(values)


def scale(polynomial, scalar):
    return trim([scalar * value for value in polynomial])


def multiply(left, right, limit=None):
    size = len(left) + len(right) - 1
    if limit is not None:
        size = min(size, limit)
    values = [0] * size
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            if left_index + right_index >= size:
                break
            values[left_index + right_index] += left_value * right_value
    return trim(values)


def power(polynomial, exponent, limit=None):
    result = [1]
    for unused in range(exponent):
        result = multiply(result, polynomial, limit)
    return result


def product(*polynomials):
    result = [1]
    for polynomial in polynomials:
        result = multiply(result, polynomial)
    return result


def derivative(polynomial):
    return trim([index * polynomial[index] for index in range(1, len(polynomial))])


def substitute(polynomial, argument, limit=None):
    result = [0]
    for coefficient in reversed(polynomial):
        result = add(multiply(result, argument, limit), [coefficient])
    return result if limit is None else result[:limit]


def coefficient(polynomial, degree):
    return polynomial[degree] if degree < len(polynomial) else 0


def recurrence_polynomials():
    return [
        scale(product([0, 1], [-1, 1], [-2, 1], [-178385058984, 83022305213], [1, 1]), 283),
        scale(product([0, 1], [-1, 1], [-2, 1],
                      [52641391340361, -45110300936020, 12196786904129]), -6),
        scale(product([-1, 1], [-2, 1],
                      [3061418504418, 202497447227007, -248299346995054, 67790996421377]), -6),
        scale(product([-2, 1], [-8342981628089220, 9123363136130610,
                               -2120668219560541, -559010679688491, 192167357944360]), -2),
        scale([-20925325629575136, 32085126156301848, -18667678380008380,
               4989857715427325, -565964885220500, 15640122296875], 5),
        scale(product([-17, 5], [-18, 5], [-21, 5], [-19, 5],
                      [-60308240636, 22509687685]), 10),
    ]


def generalized_binomial(upper, lower):
    if upper >= 0:
        return comb(upper, lower) if lower <= upper else 0
    return (-1) ** lower * comb(lower - upper - 1, lower)


def binomial_coefficients(exponent, bound, series_power=1):
    result = []
    for degree in range(bound + 1):
        numerator = 0
        for index in range(degree // 2 + 1):
            base = comb(2 * degree + series_power - 2 * index - 1, degree - 2 * index)
            base *= (-1) ** index * generalized_binomial(-exponent * (2 * degree + series_power - 2 * index), index)
            numerator += base
        quotient, remainder = divmod(series_power * numerator, degree + series_power)
        assert remainder == 0
        result.append(quotient)
    return result


def inverse_coefficients(exponent, bound):
    if exponent >= 0:
        polynomial = [0, 1, -1]
        for index in range(1, exponent + 1):
            degree = 2 * index + 1
            polynomial.extend([0] * max(0, degree + 1 - len(polynomial)))
            polynomial[degree] += (-1) ** index * comb(exponent, index)
        inverse = [0]
        for degree in range(1, bound + 2):
            residual = coefficient(substitute(polynomial, inverse, degree + 1), degree)
            inverse.append((1 if degree == 1 else 0) - residual)
        assert substitute(polynomial, inverse, bound + 2) == [0, 1]
        assert substitute(inverse, polynomial, bound + 2) == [0, 1]
    else:
        denominator = power([1, 0, -1], -exponent)
        polynomial = add([0, 1], scale(multiply([0, 0, 1], denominator), -1))
        inverse = [0]
        for degree in range(1, bound + 2):
            residual = add(substitute(polynomial, inverse, degree + 1),
                           scale(multiply([0, 1], substitute(denominator, inverse, degree + 1), degree + 1), -1))
            inverse.append(-coefficient(residual, degree))
        assert add(substitute(polynomial, inverse, bound + 2),
                   scale(multiply([0, 1], substitute(denominator, inverse, bound + 2), bound + 2), -1)) == [0]
    return inverse[1:]


def cubic_certificate():
    f_polynomial = [0, 1, -1, -1]
    h_polynomial = derivative(f_polynomial)
    numerators = [[0, 1], f_polynomial,
                  multiply(f_polynomial, add(multiply(h_polynomial, derivative(f_polynomial)),
                                             scale(multiply(derivative(h_polynomial), f_polynomial), -1)))]
    polynomials = [[0, 5, 5], [0, 11, -22], [-24, 54, -27]]
    certificate = [0]
    for lag, polynomial in enumerate(polynomials):
        shifted = substitute(polynomial, [lag - 1, 1])
        for order in range(3):
            term = product(power(f_polynomial, lag), numerators[order],
                           power(h_polynomial, 3 if order == 0 else 4 - 2 * order))
            certificate = add(certificate, scale(term, coefficient(shifted, order)))
    certificate = add(certificate, product(power(f_polynomial, 2), power(h_polynomial, 3)))
    assert certificate == [0]
    return {'result': 'zero polynomial after adding f_1^2 H_1^3',
            'operator_identity': 'L_U U_1 = -x^2; L_B B = -x',
            'recurrence_scope': 'every n >= 2 by formal Euler transport'}


def family_checks(bound):
    family_bound = min(bound, 22)
    for exponent in range(-3, 6):
        values = inverse_coefficients(exponent, family_bound)
        for series_power in range(1, 5):
            assert power(values, series_power, family_bound + 1) == binomial_coefficients(exponent, family_bound, series_power)
    cubic = inverse_coefficients(1, bound)
    for degree in range(bound + 1):
        numerator = sum(comb(degree + index, index) * comb(index, degree - index)
                        for index in range((degree + 1) // 2, degree + 1))
        assert numerator == (degree + 1) * cubic[degree]
    for degree in range(2, bound + 1):
        assert 5 * degree * (degree + 1) * cubic[degree] == 11 * degree * (2 * degree - 1) * cubic[degree - 1] + 3 * (3 * degree - 2) * (3 * degree - 4) * cubic[degree - 2]
    reciprocal = [1]
    for degree in range(1, bound + 3):
        extended = inverse_coefficients(1, degree)
        reciprocal.append(-sum(extended[degree - index] * reciprocal[index] for index in range(degree)))
    for degree in range(bound + 1):
        assert cubic[degree] == -sum(reciprocal[index] * reciprocal[degree + 2 - index]
                                    for index in range(1, degree + 3))
        if degree:
            term = Fraction(1)
            total = term
            for index in range(1, degree // 2 + 1):
                term *= Fraction((-degree + 2 * index - 2) * (-degree + 2 * index - 1), 4)
                term *= Fraction(-4, (-2 * degree + index - 1) * index)
                total += term
            assert Fraction(comb(2 * degree, degree), degree + 1) * total == cubic[degree]
    logarithmic = []
    for degree in range(1, bound + 1):
        value = sum(comb(degree + index - 1, index) * comb(index, degree - index)
                    for index in range((degree + 1) // 2, degree + 1))
        logarithmic.append(value)
        assert value == sum((index + 1) * cubic[index + 1] * reciprocal[degree - index - 1]
                            for index in range(degree))
    return {'integer_exponents': [-3, -2, -1, 0, 1, 2, 3, 4, 5],
            'series_powers': [1, 2, 3, 4], 'first_index': 0, 'last_index': family_bound,
            'cubic_alternate_sum_last_index': bound, 'cubic_recurrence_last_index': bound,
            'cubic_reciprocal_and_hypergeometric_last_index': bound,
            'A217596_initial_values': reciprocal[:10], 'A213684_initial_values': logarithmic[:10]}


def run(bound):
    polynomials = recurrence_polynomials()
    f_polynomial = [0, 1, -1, -2, 0, 1]
    h_polynomial = derivative(f_polynomial)
    numerators = [[0, 1], f_polynomial]
    for order in range(1, 5):
        step = add(multiply(h_polynomial, derivative(numerators[order])),
                   scale(multiply(derivative(h_polynomial), numerators[order]), -(2 * order - 1)))
        numerators.append(multiply(f_polynomial, step))
    w_a = [substitute(polynomial, [lag, 1]) for lag, polynomial in enumerate(polynomials)]
    w_u = [substitute(polynomial, [-1, 1]) for polynomial in w_a]
    certificate = [0]
    certificate_degree_bound = 0
    for lag in range(6):
        for order in range(6):
            term = product(power(f_polynomial, lag), numerators[order],
                           power(h_polynomial, 9 if order == 0 else 10 - 2 * order))
            certificate_degree_bound = max(certificate_degree_bound, len(term) - 1)
            certificate = add(certificate, scale(term, coefficient(w_u[lag], order)))
    assert certificate == [0]
    assert coefficient(polynomials[0], 5) == 23495312375279
    values = binomial_coefficients(2, bound)
    assert values == inverse_coefficients(2, bound)
    for degree in range(5, bound + 1):
        assert sum(substitute(polynomials[lag], [degree])[0] * values[degree - lag]
                   for lag in range(6)) == 0
    for degree in range(5):
        assert sum(substitute(polynomials[lag], [degree])[0] * values[degree - lag]
                   for lag in range(degree + 1)) == 0
    sequences = {}
    for exponent, label in [(1, 'A001002'), (2, 'A392976'), (3, 'A392977'), (-1, 'A370839')]:
        coefficients = binomial_coefficients(exponent, bound)
        assert coefficients == inverse_coefficients(exponent, bound)
        sequences[label] = coefficients[:10]
    return {
        'arithmetic': 'exact integers; polynomial coefficients in ascending degree',
        'polynomial_identity': {'result': 'zero polynomial', 'degree_bound': certificate_degree_bound,
                                'coefficients_equal_to_zero': certificate_degree_bound + 1},
        'order_five_nonzero_coefficient': 23495312375279,
        'cubic_polynomial_identity': cubic_certificate(),
        'integer_family_comparisons': family_checks(bound),
        'recurrence_checks': {'first_index': 5, 'last_index': bound, 'all_residuals': 0},
        'series_comparison': {'first_index': 0, 'last_index': bound,
                              'methods': ['binomial sums', 'recursive formal reversion'],
                              'sequences': sequences},
        'universal_scope': 'The polynomial identity and the formal-series argument establish the recurrence for every n >= 5.',
        'lean_scope': 'The pinned declaration result proves inverse identities, uniqueness, X*A=U, D-finiteness and the six-term recurrence for every n >= 5.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check the inverse-quintic polynomial identity and coefficient formulas.')
    parser.add_argument('--bound', type=int, default=40)
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    if not 5 <= arguments.bound <= 200:
        parser.error('Select a bound between 5 and 200.')
    result = run(arguments.bound)
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(result, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(text)
    print(text, end='')
