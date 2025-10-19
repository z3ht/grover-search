Grover-Search

From scratch implementation of O(√N) search algorithm developed by Lov Grover

Note: This simulation represents state as a vector of 2^n complex amplitudes, so has both time and space complexity O(2^n)

I found this project interesting because it goes against my understanding of the world. Before learning about Grover's algorithm, I thought it was impossible to search an unsorted array for a specific value in anything less than O(N).

While implementing, I developed an appreciation for the rules of quantum computing as a kind of poetry:
- Reversibility: Information cannot be discarded. I.e. `true AND false` is an illegal operation because it is impossible to infer from the resulting output `false` the original input.
- Phase Interference: The probability of a correct answer is amplified by arranging operations so that phases reinforce and cancel each other until favoring a resulting value

From this, it’s still possible to write both good and bad code.

Grover’s search is proven optimal, but with less effective interference patterns, slower solutions are also possible.

To me, this is a reminder that software is a human expression and that answers are out there waiting to be found.

Similar to classical computing, there is still a lot I don't understand about quantum computing.
