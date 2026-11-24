# Voice and humanizer revision log: report.Rmd
*2026-11-18*

Target: `analysis/report/report.Rmd` (no supplement). Passes 2
(voice, `/use_rgtvoice`) and 3 (humanizer, File mode) were
applied together, in 18 Edit calls covering every prose
paragraph outside the abstract and YAML. Nothing was staged,
committed, or rendered. The referee report for this round is
`referee-report-2026-11-18.md` in the same directory.

## Kinds of changes

1. Sentence splitting. Long sentences with subordinate
   clauses were split into two or three sentences, with a
   connective at the head of the later sentence ('However,',
   'As a result,', 'In addition,', 'That is', 'This').
   Examples: the LOCF description, the Leber dual-outcome
   paragraph, the ICH E9(R1) sentence, the Cox truth-value
   note, the reproducibility audit, the headline Results
   paragraph, and Future research items 3, 5, 6 and 7.
2. Run-on Results sentence. The headline paragraph had one
   sentence with eight inline R expressions. It now opens
   with a sentence stating the targets, then reports MMRM and
   Cox metrics in separate sentences. Every inline R
   expression is unchanged and still single-line.
3. Fragments in the Discussion. 'First, the regulatory
   requirement ... Second, ... Third, ...' were sentence
   fragments. They are now full sentences ('The first was
   ...').
4. Plain words. 'serves as' to 'is' (several), 'demonstrated'
   to 'showed', 'yielded' to 'gave', 'underpinning' to
   'basis', 'illuminate' to 'show', 'facilitate' to 'help
   with', 'ground-truth' to 'check ... against real data',
   'conducted' to 'done', 'deep regulatory precedent' to
   'long regulatory precedent', 'gold standard' to
   'standard', 'virtually every' to 'nearly every' (same
   strength), 'grossly violated' to 'badly violated'.
5. Staged openers and closers. 'This point warrants
   emphasis.' at the head of 'Note on the Cox truth value'
   was removed. It carried no information beyond the section
   heading. 'It is clear why:' became 'The reason is clear.'
   followed by the reason. 'A natural question arises: what
   is ...' became a direct question in voice ('But how
   efficient is one approach relative to the other ...?').
6. Quotation marks. Prose scare quotes "bias", "coverage"
   and "valid" became single quotes per the global rule.
   Double quotes inside backtick code spans were not touched.
7. Voice markers. 'In this paper we trace', 'We first state
   the targets', and 'We obtain the observed global CDR' put
   the collegial 'we' at section openings. Existing hedges
   ('to the best of our knowledge', 'approximately',
   'potentially', 'in general') were kept.

## Passages deliberately left alone, and why

- YAML header and abstract (hard constraint), including the
  title, which is in title case by journal convention.
- All code chunks, chunk options, and `fig.cap` strings
  (chunk options), and the `kable` caption (code).
- Display math, inline math, inline R, backtick spans,
  citation keys, and cross-references.
- ADEMP bullet lists in 'ADEMP structure', 'Design
  specifications', 'Sample-size sweep', and 'Performance
  metrics'. Their bold labels name real ADEMP elements or
  defined quantities, so they are not decorative. Only two
  over-long bullets were split.
- Bold item labels in 'Future research'. They label
  numbered items. The label of item 1 ('Additional
  covariance structures') does not match its MNAR content;
  this is referee minor comment m8 and was left for the
  remediation step.
- Section headings. All were already in sentence case; none
  were renamed, added, removed, or reordered.
- Claims flagged by the referee report were reworded at the
  same strength and not corrected. These include 'All
  participants enter with a latent global CDR of 0.5'; the
  Study 201 MMRM claim; 'nearly every pivotal AD trial';
  PMRM attributed to Jönsson et al.; 'conversion rates are
  modest' and 'heavy right-censoring'; 'real, observational
  ADNI data'; 'OS is binary'; and the placebo-dropout claim.
- Revision-history language ('as of this revision', 'now
  stops the render', 'now also requires', 'supersedes the
  2026-04-17 audit') carries information about the code
  state. It was split but kept; removing it is referee minor
  comment m1.
- The not-X-but-Y contrasts that carry information on both
  sides were kept. Examples are 'not bias against an
  independently derived ground truth', 'not pre-computed and
  pasted in', 'rather than by statistical convenience
  alone', and 'not an accident of history'.
- No sentence was merged as a duplicate. The one sentence
  deleted ('This point warrants emphasis.') is listed above.

## Verification

Script: `verify.py` in the private scratch directory. It
compares the working file with `git show HEAD:analysis/report/report.Rmd`.

| Check | Result |
|---|---|
| YAML header (including abstract) byte-identical | pass |
| Code-chunk lines, multiset (315) | identical |
| Inline R expressions, multiset (17) | identical |
| `@` citation keys, multiset (42; `\@ref` excluded) | identical |
| Numbers in prose, multiset (142) | identical |
| Math spans, multiset (45) | identical |
| Backtick spans, multiset (27) | identical |
| `\@ref` cross-references (3) | identical |
| Heading lines, in order | identical |
| Em or en dashes, ' -- ' in prose | 0 |
| Curly quotes in prose | 0 |
| British spellings | none (candidates flagged were 'exercise', 'precise', 'piecewise', 'rise', 'four', 'Otherwise') |

Sentence length. The scan counts whitespace tokens, with each
inline R, math, or code span counted as one word and
bracketed citations removed. No prose sentence exceeds 25
words. The scan still reports five hits, but all are
artifacts. Four come from bullet lists without terminal
periods, which the scanner joins into one 'sentence'. The
fifth is the ADEMP 'Aims' bullet, which counts 26 only
because the '- **Aims.**' marker is included; the sentence
itself is 24 words.

Not done: the document was not rendered, so LaTeX-level
effects of the edits were not checked. Subordinate clauses
were reduced but not eliminated. Some relative clauses
remain where splitting would have repeated a citation key or
number, which the constraints forbid.
