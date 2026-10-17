import random


def quadratic_term_count(parameter_count):
    if parameter_count < 1:
        raise ValueError("at least one parameter is required")
    return (parameter_count + 1) * (parameter_count + 2) // 2


def recommend_campaign_size(parameter_names):
    parameter_count = len(tuple(parameter_names))
    terms = quadratic_term_count(parameter_count)
    doe_trials = 3 * terms
    refinement_trials = max(12, 3 * parameter_count)
    validation_trials = max(8, 2 * parameter_count)
    return {"quadratic_terms": terms, "doe_trials": doe_trials, "refinement_trials": refinement_trials, "validation_trials": validation_trials, "max_trials": doe_trials + refinement_trials + validation_trials}


def latin_hypercube(count, bounds, seed):
    if count < 1:
        raise ValueError("count must be positive")
    rng = random.Random(seed)
    columns = {}
    for name, (low, high) in bounds.items():
        order = list(range(count))
        rng.shuffle(order)
        columns[name] = [low + (stratum + rng.random()) / count * (high - low) for stratum in order]
    return [{name: columns[name][index] for name in bounds} for index in range(count)]


def weighted_normalized_score(errors, maxima, weights):
    score = 0.0
    for name, error in errors.items():
        maximum = maxima[name]
        if maximum <= 0:
            continue
        score += weights[name] * error / maximum
    return score
