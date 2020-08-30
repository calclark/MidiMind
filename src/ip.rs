use std::{collections::HashMap, hash::Hash};

fn gen_motif_dictionary<T>(symbols: &[T]) -> HashMap<Vec<T>, u32>
where
    T: Eq + Hash + Copy,
{
    let mut dictionary = HashMap::new();
    let mut motif = vec![];

    for symbol in symbols {
        motif.push(*symbol);

        if let Some(c) = dictionary.get_mut(&motif) {
            *c += 1;
        } else {
            dictionary.insert(motif, 1);
            motif = vec![];
        }
    }

    dictionary
}

fn gen_ip_continuations<T>(motifs: HashMap<Vec<T>, u32>) -> HashMap<Vec<T>, Vec<(T, u32)>>
where
    T: Eq + Hash + Copy,
{
    let mut conts: HashMap<Vec<T>, Vec<(T, u32)>> = HashMap::new();

    for (motif, count) in motifs {
        let context = &motif[..(motif.len() - 1)];
        let symbol = &motif[motif.len() - 1];

        if let Some(c) = conts.get_mut(context) {
            c.push((*symbol, count));
        } else {
            conts.insert(context.to_vec(), vec![(*symbol, count)]);
        }
    }

    conts
}

fn normalize_continuations<T>(
    conts: HashMap<Vec<T>, Vec<(T, u32)>>,
) -> HashMap<Vec<T>, Vec<(T, f64)>>
where
    T: Eq + Hash,
{
    let mut normalized_conts = HashMap::new();

    for (context, symbols) in conts {
        let mut probabilities = vec![];
        let total_count = symbols.iter().fold(0, |acc, (_, count)| acc + count);
        for symbol in symbols {
            probabilities.push((symbol.0, symbol.1 as f64 / total_count as f64));
        }
        normalized_conts.insert(context, probabilities);
    }

    normalized_conts
}

pub fn gen_incremental_parse_tree<T>(symbols: &[T]) -> HashMap<Vec<T>, Vec<(T, f64)>>
where
    T: Eq + Hash + Copy,
{
    let motifs = gen_motif_dictionary(symbols);
    let conts = gen_ip_continuations(motifs);
    normalize_continuations(conts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn gen_ip() {
        let symbols = vec![
            "a", "b", "a", "b", "a", "b", "c", "a", "b", "d", "a", "b", "c", "d", "a", "b", "c",
            "e",
        ];
        let actual = gen_incremental_parse_tree(&symbols);
        assert_eq!(actual.len(), 4);

        for exp in vec![
            (vec![], vec![("a", 6.0 / 7.0), ("b", 1.0 / 7.0)]),
            (vec!["a"], vec![("b", 1.0)]),
            (vec!["a", "b"], vec![("c", 0.75), ("d", 0.25)]),
            (vec!["a", "b", "c"], vec![("d", 0.5), ("e", 0.5)]),
        ] {
            let conts = actual.get(&exp.0).unwrap();
            assert!(conts.len() == exp.1.len());
            for exp_cont in exp.1 {
                assert!(conts.contains(&exp_cont));
            }
        }
    }

    proptest! {
        #[test]
        fn get_ip_each_motif_continuation_probs_sum_to_one(v: Vec<i32>) {
            let probs = gen_incremental_parse_tree(&v);
            for conts in probs.values() {
                let total_prob = conts.iter().fold(0.0, |acc, cont| acc + cont.1);
                prop_assert!((total_prob - 1.0).abs() <= 0.001);
            }
        }
    }
}
