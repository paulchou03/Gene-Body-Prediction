from typing import List
import numpy as np
import pandas as pd

with open("input.fasta") as f:
    sequence = "".join(f.read().split()[1:])

states = [0, 1]
labels = {'x': 0, 'y': 1, 'z': 2, 'n': 3}

transition = np.array([
    [0.9, 0.1],
    [0.1, 0.9]
])

emission = np.array([
    [0.1, 0.1, 0.1, 0.7],
    [0.25, 0.25, 0.4, 0.1]
])

def forward_algo(seq, states, trans_prob, emiss_prob, label):
    n = len(seq)
    num_states = len(states)
    forward = np.zeros((num_states, n))
    
    for k in range(num_states):
        forward[k, 0] = (1 / num_states) * emiss_prob[k, label[seq[0]]]
    
    for i in range(1, n):
        for k in range(num_states):
            sum_prob = 0
            for l in range(num_states):
                sum_prob += forward[l, i - 1] * trans_prob[l, k]
            forward[k, i] = sum_prob * emiss_prob[k, label[seq[i]]]
        norm_factor = np.sum(forward[:, i])
        for k in range(num_states):
            forward[k, i] /= norm_factor
    
    return forward

def backward_algo(seq, states, trans_prob, emiss_prob, label):
    n = len(seq)
    num_states = len(states)
    backward = np.zeros((num_states, n))
    
    for k in range(num_states):
        backward[k, -1] = 1
    
    for i in range(n - 2, -1, -1):
        for k in range(num_states):
            sum_prob = 0
            for l in range(num_states):
                sum_prob += backward[l, i + 1] * trans_prob[k, l] * emiss_prob[l, label[seq[i + 1]]]
            backward[k, i] = sum_prob
        norm_factor = np.sum(backward[:, i])
        for k in range(num_states):
            backward[k, i] /= norm_factor
    
    return backward

def baum_welch_algo(seq, states, trans_prob, emiss_prob, label, iterations):
    num_states = len(states)
    num_symbols = len(label)
    n = len(seq)
    
    for _ in range(iterations):
        forward = forward_algo(seq, states, trans_prob, emiss_prob, label)
        backward = backward_algo(seq, states, trans_prob, emiss_prob, label)
        
        var = forward * backward
        norm_factors = np.sum(var, axis=0)
        for i in range(n):
            for k in range(num_states):
                var[k, i] /= norm_factors[i]
        
        xi = np.zeros((num_states, num_states, n - 1))
        for i in range(n - 1):
            denom = 0
            for l in range(num_states):
                for k in range(num_states):
                    denom += forward[l, i] * trans_prob[l, k] * emiss_prob[k, label[seq[i + 1]]] * backward[k, i + 1]
            for l in range(num_states):
                for k in range(num_states):
                    xi[l, k, i] = forward[l, i] * trans_prob[l, k] * emiss_prob[k, label[seq[i + 1]]] * backward[k, i + 1] / denom
        
        for l in range(num_states):
            for k in range(num_states):
                trans_prob[l, k] = np.sum(xi[l, k, :]) / np.sum(var[l, :-1])
        
        emiss_prob = np.zeros((num_states, num_symbols))
        for k in range(num_states):
            for symbol, index in label.items():
                sum = 0
                for i, s in enumerate(seq):
                    if s == symbol:
                        sum += var[k, i]
                emiss_prob[k, index] = sum
            norm_factor = np.sum(var[k, :])
            for index in range(num_symbols):
                emiss_prob[k, index] /= norm_factor
    
    return trans_prob, emiss_prob, var

transition_prob, emission_prob, var = baum_welch_algo(sequence, states, transition, emission, labels, iterations=50)

gene_body_probabilities = var[1, :]

top_50000 = sorted(
    range(len(sequence)), 
    key=lambda i: gene_body_probabilities[i], 
    reverse=True
)[:50_000]

predictions = sorted(i + 1 for i in top_50000)

with open('predictions.csv', 'w') as f:
    for interval in predictions:
        f.write(f'{interval}\n')
