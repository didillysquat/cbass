""" This script will be for the processing of any data for the 'phase 0'
publication of what I believe is the cbass_84 data set.
"""

import os
import pandas as pd

class SampleOrdinationFigure:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.input_base_dir = os.path.join(self.root_dir, 'input')
        self.distance_base_dir = os.path.join(self.input_base_dir, 'between_sample_distances')
        self.dist_path = os.path.join(self.distance_base_dir, '2019-05-13_01-00-15.556051.bray_curtis_within_clade_sample_distances.dist')
        self.pcoa_path = os.path.join(self.distance_base_dir, '2019-05-13_01-00-15.556051.PCoA_coords.csv')
        self.sample_uid_to_sample_name_dict = None
        self.sample_name_to_sample_uid_dict = {}
        self.dist_df = self._make_dist_df()

    def _make_dist_df(self):
        with open(self.dist_path, 'r') as f:
            dist_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f][1:]]

        df = pd.DataFrame(dist_data)


        self.sample_uid_to_sample_name_dict = {int(uid):name for name, uid in zip(df[0].values.tolist(), df[1].values.tolist())}
        self.sample_name_to_sample_uid_dict = {name:uid for uid, name in self.sample_uid_to_sample_name_dict.items()}
        df.drop(columns=0, inplace=True)
        df.set_index(keys=1, drop=True, inplace=True)
        df.index = df.index.astype('int')
        df.columns = df.index.values.tolist()

        return df.astype(dtype='float')

sof = SampleOrdinationFigure()
