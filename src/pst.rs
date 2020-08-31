use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
};

fn dedup_list<T>(list: &[T]) -> Vec<&T>
where
    T: Eq,
{
    let mut deduped = Vec::new();

    for item in list {
        if !deduped.contains(&item) {
            deduped.push(item);
        }
    }

    deduped
}

fn count_occurrences<T>(of: &[&T], within: &[T]) -> (u32, bool)
where
    T: Eq,
{
    let mut appears_at_end = false;
    let mut count = 0;

    for i in 0..(within.len() - of.len() + 1) {
        let mut j = 0;

        while i + j < within.len() && within[i + j] == *of[j] {
            if j == of.len() - 1 {
                count += 1;
                if i + j == within.len() - 1 {
                    appears_at_end = true;
                }
                break;
            }
            j += 1;
        }
    }

    (count, appears_at_end)
}

fn get_motif_probability<T>(motif: &[&T], within: &[T], context: Option<&[&T]>) -> f64
where
    T: Eq,
{
    let context = context.unwrap_or(&[]);

    if within.len() < motif.len() + context.len() {
        return 0.0;
    }

    let max_poss_count = if context.is_empty() {
        within.len() as u32 - motif.len() as u32 + 1
    } else {
        match count_occurrences(context, within) {
            (count, false) => count,
            (count, true) => count - 1,
        }
    };

    let mut text = context.to_vec();
    text.extend_from_slice(motif);
    let (count, _) = count_occurrences(&text, within);

    count as f64 / max_poss_count as f64
}

fn put_all_subtexts<T>(of: &[T], to: &mut HashSet<Vec<T>>)
where
    T: Eq + Hash + Copy,
{
    for i in 0..(of.len() + 1) {
        to.insert(of[i..].to_vec());
    }
}

fn get_motifs<T>(
    symbols: &[T],
    max_motif_len: usize,
    min_prob: f64,
    min_mult_fact: f64,
) -> HashSet<Vec<&T>>
where
    T: Eq + Hash,
{
    let alphabet = dedup_list(symbols);
    let mut out_motifs = HashSet::new();

    let mut poss_motifs: Vec<Vec<&T>> = alphabet
        .iter()
        .filter(|letter| get_motif_probability(&[**letter], symbols, None) >= min_prob)
        .map(|letter| vec![*letter])
        .collect();

    while !poss_motifs.is_empty() {
        let motif = poss_motifs.pop().unwrap();
        let subtext = &motif[1..];

        for letter in &alphabet {
            let cond_prob = get_motif_probability(&[*letter], symbols, Some(&motif));
            let subtext_prob = get_motif_probability(&[*letter], symbols, Some(subtext));

            let valid_mult_fact = subtext_prob == 0.0 || cond_prob / subtext_prob >= min_mult_fact;
            if cond_prob >= min_prob && valid_mult_fact {
                put_all_subtexts(&motif, &mut out_motifs);
            }
        }

        if motif.len() >= max_motif_len {
            continue;
        }

        alphabet
            .iter()
            .map(|letter| {
                let mut m = vec![*letter];
                m.extend_from_slice(&motif);
                m
            })
            .filter(|motif| get_motif_probability(&motif, symbols, None) >= min_prob)
            .for_each(|motif| poss_motifs.push(motif));
    }

    out_motifs
}

fn gen_probabilities<'a, T>(
    motifs: HashSet<Vec<&'a T>>,
    symbols: &'a [T],
    smoothing_fact: f64,
) -> HashMap<Vec<&'a T>, Vec<(&'a T, f64)>>
where
    T: Eq + Hash,
{
    let alphabet = dedup_list(symbols);
    let mut pst = HashMap::new();

    for motif in motifs {
        let mut nexts: Vec<(&T, f64)> = Vec::new();
        for letter in &alphabet {
            let prob: f64 = (1.0 - alphabet.len() as f64 * smoothing_fact)
                * get_motif_probability(&[*letter], symbols, Some(&motif))
                + smoothing_fact;
            nexts.push((letter, prob));
        }
        pst.insert(motif, nexts);
    }

    pst
}

pub fn gen_probability_suffix_tree<T>(
    symbols: &[T],
    max_motif_len: usize,
    min_prob: f64,
    min_mult_fact: f64,
    smoothing_fact: f64,
) -> HashMap<Vec<&T>, Vec<(&T, f64)>>
where
    T: Eq + Hash,
{
    let motifs = get_motifs(symbols, max_motif_len, min_prob, min_mult_fact);
    gen_probabilities(motifs, symbols, smoothing_fact)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn gen_pst() {
        let s = &["a", "b", "r", "a", "c", "a", "d", "a", "b", "r", "a"];
        let actual = gen_probability_suffix_tree(s, 2, 0.15, 2.0, 0.0);
        assert_eq!(actual.len(), 5);

        let (a, b, c, d, r) = (&"a", &"b", &"c", &"d", &"r");
        let mut maps = Vec::new();
        maps.push((
            vec![],
            vec![
                (a, 5.0 / 11.0),
                (b, 2.0 / 11.0),
                (c, 1.0 / 11.0),
                (d, 1.0 / 11.0),
                (r, 2.0 / 11.0),
            ],
        ));

        maps.push((
            vec![a],
            vec![(a, 0.0), (b, 0.5), (c, 0.25), (d, 0.25), (r, 0.0)],
        ));

        maps.push((
            vec![b],
            vec![(a, 0.0), (b, 0.0), (c, 0.0), (d, 0.0), (r, 1.0)],
        ));

        maps.push((
            vec![r],
            vec![(a, 1.0), (b, 0.0), (c, 0.0), (d, 0.0), (r, 0.0)],
        ));

        maps.push((
            vec![r, a],
            vec![(a, 0.0), (b, 0.0), (c, 1.0), (d, 0.0), (r, 0.0)],
        ));

        for (motif, exp_conts) in maps {
            let conts = actual.get(&motif).unwrap();
            assert!(conts.len() == exp_conts.len());
            for exp_cont in exp_conts {
                assert!(conts.contains(&exp_cont));
            }
        }
    }

    proptest! {
        #[test]
        fn gen_pst_each_motif_continuation_probs_sum_to_one(v: Vec<i32>) {
            let probs = gen_probability_suffix_tree(&v, 3, 0.05, 1.15, 0.02);
            for conts in probs.values() {
                let total_prob = conts.iter().fold(0.0, |acc, cont| acc + cont.1);
                prop_assert!((total_prob - 1.0).abs() <= 0.001);
            }
        }
    }
}
