// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef LDA_H
#define LDA_H

typedef struct {
    int* words;
    int* counts;
    int length;
    int total;
} document;


typedef struct {
    document* docs;
    int num_terms;
    int num_docs;
} corpus;


typedef struct {
    double alpha;
    double** log_prob_w; //beta
    int num_topics; //K
    int num_terms; //V
} lda_model;


typedef struct {
    double** class_word; // K * V
    double* class_total;
    double alpha_suffstats;
    int num_docs;
} lda_suffstats;

#endif

/*
> 1) what is the role of sufficient statistics and what does "maximum
> likelihood estimate with expected sufficient statistic" mean? Sorry if
> this sounds stupid but I am not a native speaker so this sounds sort
> of too condensed for me.
>
 - In brief, "sufficient statistics" compress the information on a
parameter (in this case, beta) to simplify its estimation (cf. Appendix
A.1 in the paper, where the sufficient statistics have a specific
meaning for each distribution). In our case, beta is a set of
multinomials (K topic-specific distributions over V terms each), and a
simple estimator for the parameters \vec\beta_k is the topic-specific
ratio between how often a term is associated to one topic versus how
often any term is associated to it (i.e., K ratios), averaged over all
documents (therefore "expected"). These "how-often" counts are not
integers, as during the E-step, only a _distribution_ over topic
associations (the variational \gamma's) for each term (in each document)
is estimated and only this can be used "for free", rather than
"decisions" on z (as for instance in a Gibbs sampler).
> what is the interpretation of \phi and \gamma variational parameters
> and how do they relate to \alpha and \beta parameters?

Regarding \alpha:

 - The estimation of _scalar_ alpha is not included in the paper. It
can, however, be derived in way similar to the one given in A.4.2 in the
LDA paper.

 - The term used as "sufficient statistic" for \alpha can be expressed
as the sum of the expected sufficient statistic E_q{\log \theta_mk} of
the document-specific topic distributions over all documents and topics.
s_\alpha = \sum_{m=1:M} \sum_{k=1:K} E_q{\log \theta_mk}

 - The expectation E_q(.) is with respect to the variational
distribution, a Dirichlet with vector parameter \vec\gamma_m. This term
becomes E_q{log \theta_mk} = digamma(\gamma_mk) - sum_k(\gamma_mk) (cf.
A.1 in the paper). So there's the connection to the rest of the algorithm.

 - In effect, the sufficient statistic "measures" the amount of
"dispersion" of parameters \theta (expressed via the variational
\gamma's) across the model. Remarkably, \alpha is only dependent on this
scalar and the two constants M and K.

Regarding \beta:

 - \sum_{m=1:M} phi_tk^(m) = E{n_kt},

 - \sum_{t=1:V} n_kt = E{n_k}

 - \beta_kt = E{n_kt} / E{n_k} \propto sum_{m=1:M} phi_kt^(m) n_mt

with E{}  the expectation over m documents and the variables n_kt
counting the occurrences of term t in topic k, n_k all occurrences of
any term in k, and n_mt the term frequency of t in document m. Note that
the phi is indexed in three dimensions while the implementation reuses a
single 2D phi array for each document (here: .^(m)).

Regarding \theta (for completeness): The \theta's can be obtained
directly from the \gamma's via \vec\theta_m = E{Dir(\vec\gamma_m)} =
\gamma_mk / \sum_k \gamma_mk.
>
> 3) in the lda-c implementation what is the meaning of lda_suffstats
> structure and class_word and class_total fields?
>
These are just the sufficient statistics for the parameter \beta:
class_word = E{n_kt}, class_total = E{n_k}

I hope this helps.

Best regards
gregor
==============================================================================
Dear Mateusz and list,

ok, I'm retrying ;)

Mateusz Berezecki wrote:
>> distributions over V terms each), and a simple estimator for the parameters
>> \vec\beta_k is the topic-specific ratio between how often a term is
>> associated to one topic versus how often any term is associated to it (i.e.,
>> K ratios), averaged over all documents (therefore "expected").
>>
>
> I am not sure if I understand this part correctly. Do you mean that
> \beta is a set of multinomials (defined in the matrix form) and that \vec\beta
> is a set of parameters (a particular distribution chosen from that matrix) ?
>
yes. \beta is a matrix of dimension K x V, and each row is a multinomial
parameter vector \vec\beta_k. Note that \phi contains the same
information but spread over terms, it's of dimension M x T_m x K, with
T_m the number of terms in document m. For completeness, I'll give some
more info on how to interpret \phi below...
> Does that ratio thing you describe refer to the following equation
> \beta_{ij} = p(w_n^j = 1 | z_n^i = 1) ?
>
>
yes. Note that my notation rather agrees with the lda-c implementation
than the LDA paper: In Eq. (9), the sum over all words of a document
(with terms w_dm^j) is equal to what I described as sum over all terms
(unique words) weighted by their frequencies, n_mt, or (n_dj if you take
the indices for documents and terms used in the paper). This shows a
lapse in my first post: Of course, the repetitions of terms in documents
also go into \phi: E{n_kt} = \sum_{m=1:M} phi_tk^(m) n_mt.
>> These
>> "how-often" counts are not integers, as during the E-step, only a
>> _distribution_ over topic associations (the variational \gamma's) for each
>> term (in each document) is estimated and only this can be used "for free",
>> rather than "decisions" on z (as for instance in a Gibbs sampler).
>>
>
> Could you please rephrase that part as I don't think I understand it.
> What do you mean that the counts are not integers?
> I think the E-step sub-sentence made the whole sentence sound
> ambiguous. Sorry if I lack understanding.
>
>
These are not really counts but rather expectations of counts, therefore
double** class_word and double* class_total, not int... By summing the
\vec\phi_t^(m) over all documents, we can get these expectations (I
didn't really mean the variational \gamma's in the first post). If you
add up (the non-integer) multinomial vectors, you get an expectation of
a total count of occurrences for each dimension of the multinomial (in
our case, topics associated with a given term t). Further explanation below.

Please consider this count interpretation a _metaphor_. I used this
because in a collapsed Gibbs sampler, there are real integer count
variables that have the same function as these float arrays. You may
have a look at Griffiths and Steyvers (PNAS, 2004) or
http://www.arbylon.net/projects/LdaGibbsSampler.java for reference.
Further, in collapsed variational Bayes (Teh et al., NIPS 2006), count
distributions (and expectations) are used explicitly.
>> - The term used as "sufficient statistic" for \alpha can be expressed as the
>> sum of the expected sufficient statistic E_q{\log \theta_mk} of the
>> document-specific topic distributions over all documents and topics.
>> s_\alpha = \sum_{m=1:M} \sum_{k=1:K} E_q{\log \theta_mk}
>>
>
> What do you mean by document-specific topic distributions over all
> documents and topics?
>
>
Consider the equation:  s_\alpha = \sum_{m=1:M} \sum_{k=1:K} E_q{\log
\theta_mk}

 - "sum" ... "over all documents and topics": \sum_{m=1:M} \sum_{k=1:K}

 - "expected suff stats": E_q{\log \theta_mk}

 - "suff stats": \log \vec\theta_m (cf. A.1 in LDA paper)

 - "expected":

 - Note that there are various quantities in the algorithm termed
sufficient statistics, some of which don't seem to be exactly the formal
sufficient statistics (e.g., cf. John Winn's thesis, 2003) but fulfil
their task of compressing the info on the parameter.

Once I'm writing anyway, another couple of remarks describing the
thoughts that sparked my understanding of LDA (and mean-field
approximations) -- maybe they help others, too:

Look into how \phi and \gamma work as priors on the z and \theta
variables in the variational approximation. They contain vague
information on their associated variables \theta and z, repectively, and
maintain a connection to global information.

For \phi, specifically: \vec\phi_t^(m) is a multinomial parameter
(distribution) of topic associations z=k to each term w=t in each
document d=m, in other words, something like p(z=k | d=m, w=t).
Averaging it over all documents (and transposing: tk -> kt) leads to the
estimate of \beta_kt = p(w=t | z=k).

It may be also worthwhile to look from where the \phi's are "fed" with
information: In the E-step, Eq. (6) or (16), respectively, global (via
the \betas) as well as local information (via the \gammas) is gathered
into \phi, while in the M-step, via Eq. (9) local information of \phi is
scattered globally.

BTW: The repetitions of terms in documents also go into \phi: E{n_kt} =
\sum_{m=1:M} phi_tk^(m) n_mt. That was missing in my first post.

Regarding \gamma:

\gamma_mk is a Dirichlet distribution over each document's \theta_m (and
represents this throughout the variational optimisation). Getting the
mean multinomial vector of this Dirichlet leads to an estimate for
\theta_mk = p(z=k | d=m).

During the E-step, the information of the document-specific topic-term
distributions is smoothly gathered into a topic distribution via Eq.
(8). This equation also reveals that some correlation exists between the
\phi and \gamma distributions of one document.

Ok, this may be redundant to what was said before, but in the hope of
better insight.

Cheers

gregor

*/
