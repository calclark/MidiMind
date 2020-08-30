use std::{collections::HashMap, hash::Hash};

fn gen_motif_dictionary<T>(symbols: &[T]) -> HashMap<Vec<&T>, u32>
where
    T: Eq + Hash,
{
    let mut dict = HashMap::new();
    let mut curr_motif = Vec::new();

    for symbol in symbols {
        curr_motif.push(symbol);

        if let Some(count) = dict.get_mut(&curr_motif) {
            *count += 1;
        } else {
            dict.insert(curr_motif, 1);
            curr_motif = Vec::new();
        }
    }

    dict
}

fn gen_continuations<T>(dict: HashMap<Vec<T>, u32>) -> HashMap<Vec<T>, Vec<(T, u32)>>
where
    T: Eq + Hash + Copy,
{
    let mut conts: HashMap<Vec<T>, Vec<(T, u32)>> = HashMap::new();

    for (motif, count) in dict {
        let context = &motif[..(motif.len() - 1)];
        let symbol = &motif[motif.len() - 1];

        if let Some(nexts) = conts.get_mut(context) {
            nexts.push((*symbol, count));
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
    let mut norm_conts = HashMap::new();

    for (context, next_symbols) in conts {
        let mut probs = Vec::new();
        let total_count = next_symbols.iter().fold(0, |acc, (_, count)| acc + count);
        for next in next_symbols {
            probs.push((next.0, next.1 as f64 / total_count as f64));
        }
        norm_conts.insert(context, probs);
    }

    norm_conts
}

pub fn gen_incremental_parse_tree<T>(symbols: &[T]) -> HashMap<Vec<&T>, Vec<(&T, f64)>>
where
    T: Eq + Hash,
{
    let dict = gen_motif_dictionary(symbols);
    let conts = gen_continuations(dict);
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
            (vec![], vec![(&"a", 6.0 / 7.0), (&"b", 1.0 / 7.0)]),
            (vec![&"a"], vec![(&"b", 1.0)]),
            (vec![&"a", &"b"], vec![(&"c", 0.75), (&"d", 0.25)]),
            (vec![&"a", &"b", &"c"], vec![(&"d", 0.5), (&"e", 0.5)]),
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
