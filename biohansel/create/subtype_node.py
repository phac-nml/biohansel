import logging
from typing import FrozenSet, List, Optional, Dict, Iterator, Type

import attr
import numpy as np
import pandas as pd
from collections import defaultdict


@attr.s
class SubtypeNode:
    """Subtype information node with links to direct ancestor and descendents if any.

    Each SubtypeNode keeps track of the genome membership and links to descendant or child subtypes and the ancestor or
    parent subtype if any.

    Attributes:
        distance: Distance level at which genome membership is defined
        genomes: Member genomes for this subtype
        children: Child/downstream subtype nodes
        parent: Parent/upstream subtype node

    """
    distance: float = attr.ib(cmp=True)
    genomes: FrozenSet[str] = attr.ib(cmp=True)
    children: List['SubtypeNode'] = attr.ib(default=attr.Factory(list), cmp=False)
    parent: Optional['SubtypeNode'] = attr.ib(default=None, repr=False, cmp=False)
    subtype: str = attr.ib(default=None, cmp=False, repr=True)
    snv_indices: pd.Int64Index = attr.ib(default=None, cmp=False)

    def build_nodes(self,
                    df_clusters: pd.DataFrame,
                    column_index: int = 0,
                    min_group_size: int = 5) -> 'SubtypeNode':
        """Build full doubly-linked hierarchical graph of SubtypeNode from a DataFrame of flat clusters a multiple levels.

        Given a full flat clusters DataFrame, the root SubtypeNode will be returned with links to all direct
        descendants (children), which can be traversed like a typical tree-like data structure.

        Note:
            This is a recursive method.

        Args:
            df_clusters: DataFrame of flat clusters at multiple levels.
            column_index: Current cluster level index to try to build SubtypeNode at.
            min_group_size: Minimum size for SubtypeNode creation.

        Returns:
            `self`: the method calling SubtypeNode instance.
        """
        if len(self.genomes) < min_group_size or column_index + 1 > df_clusters.columns.size - 1:
            return self
        distance_threshold = df_clusters.columns[column_index]
        df_clusters_subset = df_clusters.loc[self.genomes,]
        child_nodes = []
        clusters = df_clusters_subset[distance_threshold]
        cluster_counts = pd.value_counts(clusters)
        cluster_counts = cluster_counts[cluster_counts >= min_group_size]
        if cluster_counts.size == 0:
            return self

        for cluster, count in cluster_counts.iteritems():
            cluster_genomes = frozenset(clusters[clusters == cluster].index)
            if cluster_genomes == self.genomes:
                break

            child_node = SubtypeNode(distance=distance_threshold,
                                     genomes=cluster_genomes,
                                     parent=self)
            child_node.build_nodes(df_clusters=df_clusters_subset,
                                   column_index=column_index + 1,
                                   min_group_size=min_group_size)
            if child_node:
                child_nodes.append(child_node)
        if len(child_nodes) == 0:
            return self.build_nodes(df_clusters=df_clusters_subset,
                                    column_index=column_index + 1,
                                    min_group_size=min_group_size)

        self.children = child_nodes
        return self

    def drop_nodes(self, min_group_size: int = 5, depth: int = 0) -> 'SubtypeNode':
        """Drop SubtypeNodes that do not maintain an exclusive membership of a specified size.

        If an internal SubtypeNode does not meet the exclusive membership size requirement of `min_group_size`, then it
        is dropped and its children and added to its parent (i.e. grandkid nodes become child nodes). This process is
        akin to deleting items in a doubly-linked list.

        Args:
            min_group_size: Minimum group size for keeping an internal SubtypeNode.
            depth: Current depth in recursion.

        Returns:
            `self` SubtypeNode
        """
        if len(self.children) > 0:
            child_genomes = frozenset((y for x in self.children for y in x.genomes))
            count = 0
            while len(self.genomes - child_genomes) < min_group_size:
                grandkids = [y for x in self.children for y in x.children]
                for gk in grandkids:  # type: SubtypeNode
                    gk.parent = self
                self.children = grandkids
                child_genomes = frozenset((y for x in self.children for y in x.genomes))
                count += 1
            for child in self.children:
                child.drop_nodes(min_group_size=min_group_size, depth=depth + 1)

        return self

    def snvs_for_subtype(self: 'SubtypeNode', df_bin: pd.DataFrame) -> pd.Series:
        all_genomes: pd.Series = df_bin.columns
        genomes_mask = all_genomes.isin(list(self.genomes))

        row_sums_subtype = np.sum(df_bin.values[:, genomes_mask], axis=1)
        row_sums_rest = np.sum(df_bin.values[:, ~genomes_mask], axis=1)

        n_nonsubtype_genomes = all_genomes.size - len(self.genomes)
        n_subtype_genomes = len(self.genomes)
        none_subtype_all_rest_mask = (row_sums_subtype == 0) & (row_sums_rest == n_nonsubtype_genomes)
        all_subtype_none_rest_mask = (row_sums_rest == 0) & (row_sums_subtype == n_subtype_genomes)

        all_or_nothing_mask = none_subtype_all_rest_mask | all_subtype_none_rest_mask

        return df_bin.index[all_or_nothing_mask]

    def find_all_snvs(self, df_bin: pd.DataFrame, min_num_snvs: int = 2) -> 'SubtypeNode':
        snv_indices = self.snvs_for_subtype(df_bin)
        if snv_indices.size >= min_num_snvs:
            logging.info(f'Found {snv_indices.size} SNVs for differentiating SubtypeNode d={self.distance:.2e} '
                         f'N={len(self.genomes)} from everything else')
            self.snv_indices = snv_indices
        else:
            logging.info(f'Could not find enough SNVs for differentiating SubtypeNode d={self.distance:.2e} '
                         f'N={len(self.genomes)} from everything else! SNV count={snv_indices.size}')
            self.snv_indices = None

        for child in self.children:
            child.find_all_snvs(df_bin, min_num_snvs=min_num_snvs)

        return self

    def drop_snvless(self):
        dropped = []
        self._recur_drop_snvless(dropped)
        pass_count = 1
        while len(dropped) > 0:
            logging.info(f'Pass {pass_count}: Dropped {len(dropped)} SNV-less nodes')
            dropped = []
            self._recur_drop_snvless(dropped)
            pass_count += 1
        return self

    def _recur_drop_snvless(self, dropped: List['SubtypeNode'] = None):
        if dropped is None:
            dropped = []
        if self.parent and self.snv_indices is None:
            self.parent.children.remove(self)
            dropped.append(self)
            for child in self.children:
                if not child.has_children() and child.snv_indices is None:
                    continue
                child.parent = self.parent
                self.parent.children.append(child)
                child._recur_drop_snvless(dropped=dropped)
            return self
        the_children = [child for child in self.children if child.has_children() or child.snv_indices is not None]
        self.children = the_children
        for child in self.children:
            child._recur_drop_snvless(dropped=dropped)
        return self

    def has_children(self) -> bool:
        return len(self.children) > 0

    def sort_children(self) -> 'SubtypeNode':
        """Sort child nodes by number of genome from greatest to least.

        Returns:
            `self` SubtypeNode
        """
        if len(self.children) > 0:
            self.children.sort(key=lambda child: len(child.genomes), reverse=True)
            for child in self.children:
                child.sort_children()
        return self

    def assign_subtype(self,
                       depth: int = 0,
                       start_subtype: int = 1) -> 'SubtypeNode':
        """Assign biohansel subtype designations for SubtypeNode and each child SubtypeNode.

        If the parent subtype is '1' then

        - first child will have subtype '1.1'
        - second child will have subtype '1.2'
        - N-th child subtype '1.{N}'

        Notes:
            Child subtypes are incremented by one at the child level.
            Each subtype level is delimited by `.`.

        Args:
            depth: Current recursion level
            start_subtype: Starting subtype designation

        Returns:
            `self` SubtypeNode
        """
        if self.parent is None and self.snv_indices is not None:
            self.subtype = f'{start_subtype}'

        subtype = start_subtype
        for child in self.children:
            if self.subtype:
                child.subtype = f'{self.subtype}.{subtype}'
            else:
                child.subtype = f'{subtype}'
            child.assign_subtype(depth=depth + 1,
                                 start_subtype=start_subtype)
            subtype += 1
        return self

    def genome_distance_subtype(self, out: Dict[str, Dict[float, str]] = None) -> Dict[str, Dict[float, str]]:
        """From SubtypeNode graph, get 2D dict of genome to distance level to subtype designation.

        Args:
            out: Output defaultdict

        Returns:
            2D dict of genome to distance level to subtype designation.
        """
        if out is None:
            out = defaultdict(dict)
        for genome in self.genomes:
            out[genome][self.distance] = self.subtype
        for child in self.children:
            child.genome_distance_subtype(out=out)
        return out

    def node_summary(self, padright: int = 12) -> Iterator[str]:
        """Recursively get the summary of number of genomes, child genomes, distance for a node and its descendants.

        Usage:

            .. code-block:: python3

                for summary in root_subtype_node.node_summary():
                    print(summary)

            .. code-block::

                1.2          N=51 d=1.75e-03 children=132/183
                1.2.1        N=28 d=1.06e-03 children=67/95
                1.2.1.1      N=21 d=6.63e-04 children=46/67
                1.2.1.1.1    N=22 d=4.51e-04 children=24/46
                1.2.1.1.1.1  N=24 d=3.45e-04 children=0/24
                1.2.2        N=37 d=1.51e-03 children=0/37
                1.3          N=20 d=4.48e-03 children=26/46
                1.3.1        N=26 d=4.51e-04 children=0/26

        Note:
            The subtype is presented first, then the number of exclusive members in the node (`N`), then the distance
            level at which the node is defined at (`d`) and the number of `children` that the node has out of the number
            of total members for the subtype.

        Args:
            padright: Number of spaces to pad the subtype designation with so that the other information aligns nicely.

        Yields:
            Summary about node and its descendants
        """
        child_genomes = set((y for x in self.children for y in x.genomes))
        n_genomes = len(self.genomes)
        n_child_genomes = len(child_genomes)
        yield f'{(self.subtype if self.subtype else "None").ljust(padright, " ")} N={n_genomes-n_child_genomes} ' \
              f'SNVs={self.snv_indices.size if self.snv_indices is not None else None} d={self.distance:.2e} ' \
              f'children={n_child_genomes}/{n_genomes}'

        for child in self.children:
            yield from child.node_summary(padright=padright)

    @classmethod
    def from_fclusters(cls: Type['SubtypeNode'],
                       df_clusters: pd.DataFrame,
                       df_bin: pd.DataFrame,
                       min_group_size: int = 5,
                       start_subtype: int = 1,
                       min_num_snvs: int = 2) -> 'SubtypeNode':
        """Create a `SubtypeNode` tree with hierarchical subtype assignments from a table of clusters.


        Args:
            df_clusters: Table of clusters defined at multiple distance levels
            df_bin: Binary SNV matrix
            min_group_size: Minimum subtype exclusive membership size
            start_subtype: Subtype designation to start at. If `0`, then the first subtype will be `0` and the next `1`, `2`, and so on.
            min_num_snvs: Minimum number of SNVs for a subtype

        Returns:
            Root `SubtypeNode` with links to descendants if any and hierarchical subtypes assigned.
        """

        root_genomes = frozenset(df_clusters.index)
        root_distance = df_clusters.columns[0]

        obj = cls(distance=root_distance,
                  genomes=root_genomes,
                  parent=None)
        logging.info(f'Initialized root SubtypeNode with distance={obj.distance} and N={len(obj.genomes)} members.')
        logging.info(f'Building full tree of SubtypeNodes from clusters table of size {df_clusters.shape} and minimum '
                     f'group size of {min_group_size}')
        obj.build_nodes(df_clusters=df_clusters, min_group_size=min_group_size)
        logging.info(f'Dropping internal SubtypeNodes that do not meet minimum size requirement of {min_group_size}')
        obj.drop_nodes(min_group_size=min_group_size)

        row_sums = np.sum(df_bin.values, axis=1)
        df_bin_filtered = df_bin[row_sums >= min_group_size]

        logging.info(f'Filtered binary SNV matrix from N={df_bin.shape[0]} rows to N={df_bin_filtered.shape[0]} rows '
                     f'based on minimum group size of {min_group_size}. Only SNVs present in at least {min_group_size} '
                     f'genomes will be used to find SNVs differentiating subgroups.')
        logging.info(f'Finding SNVs that differentiate genomes belonging to each SubtypeNode from binary SNV matrix '
                     f'with size {df_bin_filtered.shape}')
        obj.find_all_snvs(df_bin, min_num_snvs=min_num_snvs)
        logging.info(f'Dropping SubtypeNodes that do not contain enough SNVs')
        obj.drop_snvless()
        logging.info(f'Sorting all children for all SubtypeNodes by size from largest to smallest')
        obj.sort_children()
        logging.info(f'Assigning hierarchical subtype designations for all SubtypeNodes')
        obj.assign_subtype(start_subtype=start_subtype)

        return obj

    def leaves(self) -> Iterator['SubtypeNode']:
        """Get iterator for all leaf SubtypeNodes"""
        if len(self.children) == 0:
            yield self
        for c in self.children:
            yield from c.leaves()

    def root(self) -> 'SubtypeNode':
        """Get root SubtypeNode"""
        if self.parent:
            return self.parent.root()
        else:
            return self

    def _summary_table_row(self):
        child_genomes = set((y for x in self.children for y in x.genomes))
        n_genomes = len(self.genomes)
        n_child_genomes = len(child_genomes)
        return dict(subtype=self.subtype,
                    exclusive_genomes=n_genomes - n_child_genomes,
                    child_genomes=n_child_genomes,
                    total_genomes=n_genomes,
                    n_snvs=(self.snv_indices.size if self.snv_indices is not None else 0),
                    distance_level=self.distance)

    def _recur_get_summary_table_rows(self, rows):
        rows.append(self._summary_table_row())
        for child in self.children:
            child._recur_get_summary_table_rows(rows=rows)

    def summary_table(self):
        table_rows = []
        self._recur_get_summary_table_rows(table_rows)
        df = pd.DataFrame(table_rows)

        cols = [x.strip() for x in '''
        subtype
        exclusive_genomes
        child_genomes
        total_genomes
        n_snvs
        distance_level
        '''.strip().split('\n')]

        df = df[cols]
        return df

    def subtype_level_table(self, index=None, view_in_phandango=True):
        gds = self.genome_distance_subtype()
        genome_levels_subtype = {g: {i: ds[d] for i, d in enumerate(sorted(list(ds.keys()), reverse=True))}
                                 for g, ds in gds.items()}
        df = pd.DataFrame(genome_levels_subtype).transpose()
        if pd.isnull(df[0]).all():
            df.drop(columns=[0], inplace=True)

        if index is not None:
            df = df.loc[index,]

        if view_in_phandango:
            df.columns = [f'L{c}:o' for c in df.columns]

        return df
