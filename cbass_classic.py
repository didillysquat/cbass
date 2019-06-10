""" This script will be for the processing of any data for the 'phase 0'
publication of what I believe is the cbass_classic data set.
"""

import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
# NB the pip cartopy install seems to be broken as it doesn't install the required libararies.
# The solution was to install using conda. conda install cartopy.
# I then had to downgrade shapely to 1.5.17. pip install shapely==1.5.17
from cartopy.mpl.gridliner import Gridliner
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy
import matplotlib.gridspec as gridspec
from matplotlib import collections, patches
import sys
import random
from matplotlib.patches import Rectangle, Polygon, Arrow
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import numpy as np

class SampleOrdinationFigure:
    def __init__(self, plot_seqs_by='type', distance_plotting_format='only_samples', seq_type_arrangement='by_site'):
        self.plot_seqs_by = plot_seqs_by
        self.seq_type_arrangement = seq_type_arrangement
        self.distance_plotting_format = distance_plotting_format
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.input_base_dir = os.path.join(self.root_dir, 'classic', 'input')

        # sequences abundances
        self.seq_abund_relative_path = os.path.join(self.input_base_dir, '52_DBV_21022019_2019-06-09_00-56-14.628598.seqs.relative.txt')
        self.seq_rel_df = self._make_seq_rel_df()
        self.ordered_seq_names = self._get_ordered_seq_names()
        self.seq_pal = self._get_colour_list()
        self.seq_color_dict = self._get_seq_color_dict()

        # profile abundances
        self.profile_abund_relative_path = os.path.join(self.input_base_dir, '52_DBV_21022019_2019-06-09_00-56-14.628598.profiles.relative.txt')
        self.profile_uid_to_profile_name_dict = {}
        self.profile_name_to_profile_uid_dict = {}
        self.prof_rel_df = self._make_prof_rel_df()
        self.ordered_prof_names = self._get_ordered_prof_names()
        self.prof_pal = self._get_prof_pal()
        self.prof_color_dict = self._get_prof_color_dict()

        # Btwn sample distances
        self.distance_base_dir_samples = os.path.join(self.input_base_dir, 'between_sample_distances')
        self.dist_path_samples = os.path.join(self.distance_base_dir_samples, 'A','2019-06-09_00-56-14.628598.bray_curtis_sample_distances_A.dist')
        self.pcoa_samples_path = os.path.join(self.distance_base_dir_samples, 'A','2019-06-09_00-56-14.628598.bray_curtis_samples_PCoA_coords_A.csv')
        self.sample_uid_to_sample_name_dict = None
        self.sample_name_to_sample_uid_dict = None
        self.sample_dist_df = self._make_sample_dist_df()
        self.sample_pcoa_df = self._make_sample_pcoa_df()

        # Btwn type distances
        self.distance_base_dir_types = os.path.join(self.input_base_dir, 'between_profile_distances', 'A')
        self.dist_path_types = os.path.join(self.distance_base_dir_types,
                                              '2019-06-09_00-56-14.628598.bray_curtis_within_clade_profile_distances_A.dist')
        self.pcoa_types_path = os.path.join(self.distance_base_dir_types, '2019-06-09_00-56-14.628598.bray_curtis_profiles_PCoA_coords_A.csv')
        self.type_uid_to_type_name_dict = None
        self.type_name_to_type_uid_dict = None
        self.type_dist_df = self._make_type_dist_df()
        self.type_pcoa_df = self._make_type_pcoa_df()

        # metainfo
        self.meta_path = os.path.join(self.root_dir, '84','input', 'meta_info.xlsx')
        self.meta_df = self._make_meta_df()

        self.fig_out_path = os.path.join(self.root_dir, 'classic', 'figures')
        os.makedirs(self.fig_out_path, exist_ok=True)
        # self.gis_input_base_path = os.path.join(self.input_base_dir, 'gis')

        # Transcriptomics
        # self.transcript_input_base_path = os.path.join(self.input_base_dir, 'transcripts')
        # self.transcript_kallisto_matrix_path = os.path.join(self.transcript_input_base_path, 'spis_kallisto_TPM_norm_matrix.txt')
        # self.transcript_output_path = os.path.join(self.transcript_input_base_path, 'outputs')
        # self.host_kallisto_df = self._make_kallisto_df()

        # Figure setup
        # have it set up so that we can transfer this over to the cbass_84 script.
        self.fig = plt.figure(figsize=(8, 8))
        if self.distance_plotting_format == 'only_samples':
            self.gs = gridspec.GridSpec(2, 2)
            self.large_map_ax = plt.subplot(self.gs[:1, :1], projection=ccrs.PlateCarree(), zorder=1)
            self.pc1_pc2_sample_dist_ax = plt.subplot(self.gs[1:2, :1])
            self.pc1_pc3_sample_dist_ax = plt.subplot(self.gs[1:2, 1:2])
            self.sub_plot_gs = gridspec.GridSpecFromSubplotSpec(2, 55, subplot_spec=self.gs[0:1, 1:2])
            self.sub_plot_seqs_axarr = [plt.subplot(self.sub_plot_gs[0:1, 0:55]),
                                        plt.subplot(self.sub_plot_gs[1:2, 0:29]),
                                        plt.subplot(self.sub_plot_gs[1:2, 29:])]
        elif self.distance_plotting_format == 'sample_type':
            relative_width_large_map = 16/38
            relative_width_sub_plot_gs = 1/2
            gs_rows = 4
            # 8 for each plot, 2 for gap in between with 3 gaps = 8*4 + 3*2 = 38
            plot_width = 8
            gap_width = 2
            num_dist_plots = 4
            gs_cols = (num_dist_plots*plot_width) + ((num_dist_plots-1)*gap_width)
            self.gs = gridspec.GridSpec(gs_rows, gs_cols)

            # self.large_map_ax = plt.subplot(self.gs[:3, :int(relative_width_large_map*gs_cols)], projection=ccrs.PlateCarree(), zorder=1)
            self.large_map_ax = plt.subplot(self.gs[:3, :int(relative_width_large_map*gs_cols)], projection=ccrs.PlateCarree(), zorder=1)

            self.pc1_pc2_sample_dist_ax = plt.subplot(self.gs[3:4, :plot_width])
            self.pc1_pc3_sample_dist_ax = plt.subplot(self.gs[3:4, plot_width + gap_width:(plot_width * 2) + gap_width])
            self.pc1_pc2_type_dist_ax = plt.subplot(self.gs[3:4, (plot_width * 2) + (gap_width * 2):(plot_width * 3) + (gap_width * 2)])
            self.pc1_pc3_type_dist_ax = plt.subplot(self.gs[3:4, (plot_width * 3) + (gap_width * 3):])
            self.sub_plot_gs = gridspec.GridSpecFromSubplotSpec(2, 55, subplot_spec=self.gs[0:3, int(gs_cols*relative_width_sub_plot_gs):gs_cols])
            self.sub_plot_seqs_axarr = [plt.subplot(self.sub_plot_gs[0:1, 0:55]),
                                        plt.subplot(self.sub_plot_gs[1:2, 0:29]),
                                        plt.subplot(self.sub_plot_gs[1:2, 29:])]
        elif self.distance_plotting_format == 'sample_type_by_site':
            relative_width_large_map = 16 / 38
            relative_width_sub_plot_gs = 1 / 2
            gs_rows = 4
            # 8 for each plot, 2 for gap in between with 3 gaps = 8*4 + 3*2 = 38
            plot_width = 8
            gap_width = 2
            num_dist_plots = 4
            gs_cols = (num_dist_plots * plot_width) + ((num_dist_plots - 1) * gap_width)
            self.gs = gridspec.GridSpec(gs_rows, gs_cols)

            # self.large_map_ax = plt.subplot(self.gs[:3, :int(relative_width_large_map*gs_cols)], projection=ccrs.PlateCarree(), zorder=1)
            self.large_map_ax = plt.subplot(self.gs[:3, :int(relative_width_large_map * gs_cols)],
                                            projection=ccrs.PlateCarree(), zorder=1)

            self.pc1_pc2_sample_dist_ax = plt.subplot(self.gs[3:4, :plot_width])
            self.pc1_pc3_sample_dist_ax = plt.subplot(self.gs[3:4, plot_width + gap_width:(plot_width * 2) + gap_width])
            self.pc1_pc2_type_dist_ax = plt.subplot(
                self.gs[3:4, (plot_width * 2) + (gap_width * 2):(plot_width * 3) + (gap_width * 2)])
            self.pc1_pc3_type_dist_ax = plt.subplot(self.gs[3:4, (plot_width * 3) + (gap_width * 3):])
            site_plot_width = 21
            site_gap_width = 4
            sub_plot_gs_width = ((2*site_plot_width) + site_gap_width)
            self.sub_plot_gs = gridspec.GridSpecFromSubplotSpec(6, sub_plot_gs_width, subplot_spec=self.gs[0:3, int(
                gs_cols * relative_width_sub_plot_gs):gs_cols])
            self.sub_plot_seqs_axarr = [plt.subplot(self.sub_plot_gs[0:2, :site_plot_width]),
                                        plt.subplot(self.sub_plot_gs[0:2, site_plot_width + site_gap_width:]),
                                        plt.subplot(self.sub_plot_gs[2:4, :site_plot_width]),
                                        plt.subplot(self.sub_plot_gs[2:4, site_plot_width + site_gap_width:]),
                                        plt.subplot(self.sub_plot_gs[4:5, :]),
                                        plt.subplot(self.sub_plot_gs[5:6, :]) ]

        # self.zoom_map_ax = plt.subplot(self.gs[:1, 1:2], projection=ccrs.PlateCarree(), zorder=1)

        self.sites = ['eilat', 'kaust', 'exposed', 'protected']
        self.sites_location_dict = {'eilat': (34.934402, 29.514673), 'kaust': (38.960234, 22.303411), 'exposed': (39.04470, 22.270300), 'protected':(39.051982,22.26900)}
        self.site_color_dict = {'eilat':'white', 'kaust': 'white', 'exposed': 'white', 'protected':'black'}

        # self.sub_plot_profiles_axarr = [plt.subplot(self.sub_plot_gs[1:2, 0:1]), plt.subplot(self.sub_plot_gs[3:4, 0:1])]
        self.site_marker_dict = {'eilat': 'o', 'kaust': '^', 'protected': '+', 'exposed': 's' }

    def _get_prof_pal(self):
        return ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                         self.create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000, num_cols=6,
                                                 time_out_iterations=100000)]

    def _get_prof_color_dict(self):
        temp_col_dict = {}
        for i, prof_name in enumerate(self.ordered_prof_names):
            temp_col_dict[prof_name] = self.prof_pal[i]
        return temp_col_dict

    def _get_ordered_prof_names(self):
        return self.prof_rel_df.sum().sort_values(ascending=False).index.values.tolist()

    def _get_seq_color_dict(self):
        temp_col_dict = {}
        for i, seq in enumerate(self.ordered_seq_names):
            temp_col_dict[seq] = self.seq_pal[i]
        return temp_col_dict

    def _get_ordered_seq_names(self):
        return self.seq_rel_df.sum().sort_values(ascending=False).index.values.tolist()

    def _make_prof_rel_df(self):
        with open(self.profile_abund_relative_path, 'r') as f:
            prof_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f]]

        df = pd.DataFrame(prof_data)
        self.profile_uid_to_profile_name_dict = {uid:name for uid, name in zip(df.iloc[0,2:].values.tolist(), df.iloc[6, 2:].values.tolist())}
        self.profile_name_to_profile_uid_dict = {name:uid for uid, name in self.profile_uid_to_profile_name_dict.items()}
        df.drop(index=list(range(6)), inplace=True)
        df.drop(columns=1, inplace=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:-6,:]
        df.set_index('ITS2 type profile', drop=True, inplace=True)
        df.index = df.index.astype('int')
        return df.astype('float')

    def _make_seq_rel_df(self):
        with open(self.seq_abund_relative_path, 'r') as f:
            seq_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f]]

        df = pd.DataFrame(seq_data)
        df.iat[0,0] = 'sample_uid'
        df.columns = df.iloc[0]
        df.drop(index=0, inplace=True)
        df.drop(columns='sample_name', inplace=True)
        df.set_index('sample_uid', drop=True, inplace=True)
        df = df[:-3]
        df.index = df.index.astype('int')
        # Get rid of all of the superflous columns only leaving the seq rel counts
        df = df.iloc[:,20:]

        return df.astype('float')

    def _make_sample_dist_df(self):
        with open(self.dist_path_samples, 'r') as f:
            dist_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f][1:]]

        df = pd.DataFrame(dist_data)
        self.sample_uid_to_sample_name_dict = {int(uid):name for name, uid in zip(df[0].values.tolist(), df[1].values.tolist())}
        self.sample_name_to_sample_uid_dict = {name:uid for uid, name in self.sample_uid_to_sample_name_dict.items()}
        df.drop(columns=0, inplace=True)
        df.set_index(keys=1, drop=True, inplace=True)
        df.index = df.index.astype('int')
        df.columns = df.index.values.tolist()

        return df.astype(dtype='float')

    def _make_sample_pcoa_df(self):
        df = pd.read_csv(filepath_or_buffer=self.pcoa_samples_path, sep=',', header=0, index_col=1, skipfooter=1)
        df.drop(columns='sample', inplace=True)
        df.index = df.index.astype('int')
        return df

    def _make_type_dist_df(self):
        with open(self.dist_path_types, 'r') as f:
            dist_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f][1:]]

        df = pd.DataFrame(dist_data)
        self.type_uid_to_type_name_dict = {int(uid): name for name, uid in
                                               zip(df[0].values.tolist(), df[1].values.tolist())}
        self.type_name_to_type_uid_dict = {name: uid for uid, name in self.type_uid_to_type_name_dict.items()}
        df.drop(columns=0, inplace=True)
        df.set_index(keys=1, drop=True, inplace=True)
        df.index = df.index.astype('int')
        df.columns = df.index.values.tolist()

        return df.astype(dtype='float')

    def _make_type_pcoa_df(self):
        df = pd.read_csv(filepath_or_buffer=self.pcoa_types_path, sep=',', header=0, index_col=1, skipfooter=1)
        df.drop(columns='sample', inplace=True)
        df.index = df.index.astype('int')
        return df

    def _make_meta_df(self):
        df = pd.read_excel(self.meta_path, sheet_name=0, header=1, index_col=False, drop_row=0 )
        df.set_index('sample_name', drop=True, inplace=True)
        # make sure that each of the sample names in the distance file are found in the meta file
        sample_names_in_meta = df.index.values.tolist()
        ordered_sample_tups_list = []
        for sample_name, sample_uid in self.sample_name_to_sample_uid_dict.items():
            if sample_name not in sample_names_in_meta:
                print(f'Sample name: {sample_name} not found')
            else:
                ordered_sample_tups_list.append((sample_name, sample_uid))
        wanted_sample_names = [tup[0] for tup in ordered_sample_tups_list]
        # fix the bad coordinate system and then convert the coordinate to float for theses columns
        for i, sample_name in enumerate(df.index.values.tolist()): # for each row
            current_val_lat = df['collection_latitude'][sample_name]
            if not isinstance(current_val_lat, float):
                if 'N' in current_val_lat:
                    df.iat[i, 9] = float(current_val_lat[2:-1])
                    current_val_lon = df['collection_longitude'][sample_name]
                    df.iat[i, 10] = float(current_val_lon[2:-1])
        df = df.loc[wanted_sample_names,:]
        # make a new columns that is site name
        site_name = []
        for i, sample_name in enumerate(df.index.values.tolist()):
            if df['collection_latitude'][sample_name] == 29.514673:
                site_name.append('eilat')
            elif df['collection_latitude'][sample_name] == 22.26302:
                site_name.append('protected')
            elif df['collection_latitude'][sample_name] == 22.26189:
                site_name.append('exposed')
            elif df['collection_latitude'][sample_name] == 22.303411:
                site_name.append('kaust')
            else:
                sys.exit('Mismatch in latitude')
        df['site'] = site_name
        # make a new column that is classic or cbass
        cbass_classic_list = []
        for i, sample_name in enumerate(df.index.values.tolist()):
            if 'ST' in sample_name:
                cbass_classic_list.append('CBASS')
            elif 'LT' in sample_name:
                cbass_classic_list.append('Classic')
            else:
                cbass_classic_list.append('n/a')
        df['cbass_classic'] = cbass_classic_list

        return df

    @staticmethod
    def _get_colour_list():
        return [
            "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900",
            "#0000A6",
            "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400",
            "#4FC601",
            "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
            "#B903AA",
            "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101",
            "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
            "#0CBD66",
            "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
            "#BEC459",
            "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
            "#FF913F",
            "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7",
            "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
            "#201625",
            "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
            "#CB7E98",
            "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
            "#806C66",
            "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5",
            "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
            "#353339",
            "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
            "#001325",
            "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
            "#8D8546",
            "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F",
            "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
            "#A45B02",
            "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
            "#F4D749",
            "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
            "#C6DC99",
            "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A",
            "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
            "#A38469",
            "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
            "#E4FFFC",
            "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
            "#A76F42",
            "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2",
            "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
            "#868E7E",
            "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
            "#00B57F",
            "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]

    def create_colour_list(
            self, sq_dist_cutoff=None, mix_col=None, num_cols=50,
            time_out_iterations=10000, avoid_black_and_white=True):
        new_colours = []
        min_dist = []
        attempt = 0
        while len(new_colours) < num_cols:
            attempt += 1
            # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
            if attempt > time_out_iterations:
                sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                         'and have not been able to find a colour that fits into your defined colour space.\n'
                         'Please lower the number of colours you are trying to find, '
                         'the minimum distance between them, or both.'.format(attempt))
            if mix_col:
                r = int((random.randint(0, 255) + mix_col[0]) / 2)
                g = int((random.randint(0, 255) + mix_col[1]) / 2)
                b = int((random.randint(0, 255) + mix_col[2]) / 2)
            else:
                r = random.randint(0, 255)
                g = random.randint(0, 255)
                b = random.randint(0, 255)

            # now check to see whether the new colour is within a given distance
            # if the avoids are true also
            good_dist = True
            if sq_dist_cutoff:
                dist_list = []
                for i in range(len(new_colours)):
                    distance = (new_colours[i][0] - r) ** 2 + (new_colours[i][1] - g) ** 2 + (
                                new_colours[i][2] - b) ** 2
                    dist_list.append(distance)
                    if distance < sq_dist_cutoff:
                        good_dist = False
                        break
                # now check against black and white
                d_to_black = (r - 0) ** 2 + (g - 0) ** 2 + (b - 0) ** 2
                d_to_white = (r - 255) ** 2 + (g - 255) ** 2 + (b - 255) ** 2
                if avoid_black_and_white:
                    if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                        good_dist = False
                if dist_list:
                    min_dist.append(min(dist_list))
            if good_dist:
                new_colours.append((r, g, b))
                attempt = 0

        return new_colours

sof = SampleOrdinationFigure(distance_plotting_format='sample_type_by_site')
