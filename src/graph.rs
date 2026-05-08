use color_eyre::Result;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(usize);

impl NodeId {
    pub fn get_id(&self) -> usize {
        self.0
    }
}

#[derive(Clone, Debug)]
pub struct Graph {
    node_num: usize,
    adjacency: HashMap<NodeId, HashSet<NodeId>>,
}

impl Graph {
    fn new() -> Self {
        Graph {
            node_num: 0,
            adjacency: HashMap::new(),
        }
    }

    fn new_node(&mut self) -> NodeId {
        let node_id = NodeId(self.node_num);
        self.node_num += 1;
        node_id
    }

    fn add_edge(&mut self, u: NodeId, v: NodeId) -> Result<()> {
        if u == v {
            return Err(color_eyre::eyre::eyre!("自己ループは許可されていません"));
        }
        self.adjacency
            .entry(u)
            .or_default()
            .insert(v);
        self.adjacency
            .entry(v)
            .or_default()
            .insert(u);
        Ok(())
    }

    pub fn neighbors(&self, node: NodeId) -> Option<&HashSet<NodeId>> {
        self.adjacency.get(&node)
    }

    pub fn iter_nodes(&self) -> impl Iterator<Item = NodeId> + '_ {
        (0..self.node_num).map(NodeId)
    }
    pub fn get_node_num(&self) -> usize {
        self.node_num
    }
}

pub fn generate_graph(n: usize, edges: &[(usize, usize)]) -> Result<Graph> {
    let mut graph = Graph::new();
    let node_ids = {
        let mut ids = Vec::new();
        for _ in 0..n {
            ids.push(graph.new_node());
        }
        ids
    };
    for &(u, v) in edges {
        let (Some(u_id), Some(v_id)) = (node_ids.get(u), node_ids.get(v)) else {
            return Err(color_eyre::eyre::eyre!(
                "ノードID {} または {} はグラフのノード数 {} を超えています",
                u,
                v,
                n
            ));
        };
        graph.add_edge(*u_id, *v_id)?;
    }
    Ok(graph)
}

/// 連結成分分解
/// グラフが連結でない場合、複数の連結成分に分割する
fn connected_components(graph: &Graph) -> Result<Vec<Graph>> {
    let mut visited = HashSet::new();
    let mut components = Vec::new();

    for node in graph.iter_nodes() {
        if !visited.contains(&node) {
            let mut stack = vec![node];
            let mut nodes_in_component = Vec::new();
            while let Some(current) = stack.pop() {
                if visited.insert(current) {
                    nodes_in_component.push(current);
                    if let Some(neighbors) = graph.neighbors(current) {
                        for neighbor in neighbors {
                            stack.push(*neighbor);
                        }
                    }
                }
            }
            let mut component = Graph::new();
            let mut node_mapping = HashMap::new();
            let sorted_nodes_in_component = {
                let mut sorted = nodes_in_component.clone();
                sorted.sort_by_key(|n| n.get_id());
                sorted
            };
            for node in sorted_nodes_in_component {
                let new_id = component.new_node();
                node_mapping.insert(node, new_id);
            }
            // node_mappingをVectorに変換して、元のノードIDの順序を保つ
            let sorteld_node_mapping = {
                let mut sorted = node_mapping.iter().collect::<Vec<_>>();
                sorted.sort_by_key(|(node, _)| node.get_id());
                sorted
            };
            for (node, &new_id) in sorteld_node_mapping {
                if let Some(neighbors) = graph.neighbors(*node) {
                    for &neighbor in neighbors {
                        if let Some(&new_neighbor_id) = node_mapping.get(&neighbor) {
                            component.add_edge(new_id, new_neighbor_id)?;
                        }
                    }
                }
            }
            components.push(component);
        }
    }
    Ok(components)
}

/// 次数が1のノードを削除する
fn remove_degree_one_nodes(graph: &Graph) -> Result<Graph> {
    let nodes = graph
        .iter_nodes()
        .filter(|node| {
            if let Some(neighbors) = graph.neighbors(*node) {
                neighbors.len() > 1
            } else {
                false
            }
        })
        .collect::<Vec<_>>();
    // 新しいグラフを作成
    let mut new_graph = Graph::new();
    // nodeからnew_nodeへのマッピングを作成
    let node_to_newnode = nodes
        .iter()
        .map(|node| (node, new_graph.new_node()))
        .collect::<HashMap<_, _>>();
    // 新しいグラフに辺を追加
    for node in nodes.iter() {
        if let Some(neighbors) = graph.neighbors(*node)
            && let Some(&new_node) = node_to_newnode.get(node)
        {
            for neighbor in neighbors {
                if let Some(&new_neighbor) = node_to_newnode.get(&neighbor) {
                    new_graph.add_edge(new_node, new_neighbor)?;
                }
            }
        }
    }
    Ok(new_graph)
}

/// 次数が1のノードを削除する処理を繰り返す
fn remove_degree_one_nodes_iteratively(graph: &Graph) -> Result<Graph> {
    let mut current_graph = graph.clone();
    loop {
        let new_graph = remove_degree_one_nodes(&current_graph)?;
        if new_graph.node_num == current_graph.node_num {
            break Ok(new_graph);
        }
        current_graph = new_graph;
    }
}

/// 連結成分分解をしたうえで、次数が1のノードを削除する処理を繰り返す
pub fn simplify_graph(graph: &Graph) -> Result<Vec<Graph>> {
    let components = connected_components(graph)?;
    components
        .iter()
        .map(remove_degree_one_nodes_iteratively)
        .collect()
}

impl std::fmt::Display for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Node count: {}", self.node_num)?;
        for node in self.iter_nodes() {
            if let Some(neighbors) = self.neighbors(node) {
                let sorted_neighbors = {
                    let mut sorted = neighbors.iter().cloned().collect::<Vec<_>>();
                    sorted.sort_by_key(|n| n.get_id());
                    sorted
                };
                write!(f, "{}: ", node)?;
                for neighbor in sorted_neighbors {
                    write!(f, "{} ", neighbor)?;
                }
                writeln!(f)?;
            } else {
                writeln!(f, "{}: No neighbors", node)?;
            }
        }
        Ok(())
    }
}

// test
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_connected_components() -> Result<()> {
        let graph = generate_graph(6, &[(0, 4), (1, 3), (3, 2), (2, 1)])?;
        let components = connected_components(&graph)?;
        assert_eq!(components.len(), 3);
        let components_size: Vec<usize> = components.iter().map(|c| c.node_num).collect();
        assert!(components_size.contains(&1));
        assert!(components_size.contains(&2));
        assert!(components_size.contains(&3));
        Ok(())
    }
}
