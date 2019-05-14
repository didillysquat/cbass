""" This script will be for the processing of any data for the 'phase 0'
publication of what I believe is the cbass_84 data set.
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
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import numpy as np

class SampleOrdinationFigure:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.input_base_dir = os.path.join(self.root_dir, 'input')
        self.distance_base_dir = os.path.join(self.input_base_dir, 'between_sample_distances')
        self.dist_path = os.path.join(self.distance_base_dir, '2019-05-13_08-54-39.673986.bray_curtis_within_clade_sample_distances.dist')
        self.pcoa_path = os.path.join(self.distance_base_dir, '2019-05-13_08-54-39.673986.PCoA_coords.csv')
        self.meta_path = os.path.join(self.input_base_dir, 'meta_info.xlsx')
        self.profile_abund_relative_path = os.path.join(self.input_base_dir, '52_DBV_21022019_2019-04-21_09-11-11.379408.profiles.relative.txt')
        self.seq_abund_relative_path = os.path.join(self.input_base_dir, '52_DBV_21022019_2019-04-21_09-11-11.379408.seqs.relative.txt')
        self.sample_uid_to_sample_name_dict = None
        self.sample_name_to_sample_uid_dict = {}
        self.dist_df = self._make_dist_df()
        self.pcoa_df = self._make_pcoa_df()
        self.meta_df = self._make_meta_df()
        self.seq_rel_df = self._make_seq_rel_df()
        self.ordered_seq_names = self._get_ordered_seq_names()
        self.seq_pal = self._get_colour_list()
        self.seq_color_dict = self._get_seq_color_dict()
        self.profile_uid_to_profile_name_dict = {}
        self.profile_name_to_profile_uid_dict = {}
        self.prof_rel_df = self._make_prof_rel_df()
        self.ordered_prof_uids = self._get_ordered_prof_uids()
        self.prof_pal = self._get_prof_pal()
        self.prof_color_dict = self._get_prof_color_dict()
        self.fig = plt.figure(figsize=(8, 12))
        self.gs = gridspec.GridSpec(2, 2)
        self.large_map_ax = plt.subplot(self.gs[:1, :1], projection=ccrs.PlateCarree(), zorder=1)
        # self.zoom_map_ax = plt.subplot(self.gs[:1, 1:2], projection=ccrs.PlateCarree(), zorder=1)
        self.pc1_pc2_ax = plt.subplot(self.gs[1:2, :1])
        self.sites = ['eilat', 'kaust', 'exposed', 'protected']
        self.sites_location_dict = {'eilat': (34.934402, 29.514673), 'kaust': (38.960234, 22.303411), 'exposed': (39.04878, 22.26189), 'protected':(39.05165,22.26302)}
        self.site_color_dict = {'eilat':'white', 'kaust': 'white', 'exposed': 'white', 'protected':'black'}
        self.pc1_pc3_ax = plt.subplot(self.gs[1:2, 1:2])
        self.sub_plot_gs = gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec=self.gs[0:1, 1:2])
        self.sub_plot_seqs_axarr = [plt.subplot(self.sub_plot_gs[0:1,0:1]), plt.subplot(self.sub_plot_gs[1:2,0:1])]
        # self.sub_plot_profiles_axarr = [plt.subplot(self.sub_plot_gs[1:2, 0:1]), plt.subplot(self.sub_plot_gs[3:4, 0:1])]
        self.site_marker_dict = {'eilat': 'o', 'kaust': '^', 'protected': '+', 'exposed': 's' }

    def _get_prof_pal(self):
        return ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                         self.create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000, num_cols=50,
                                                 time_out_iterations=10000)]

    def _get_prof_color_dict(self):
        temp_col_dict = {}
        for i, prof_uid in enumerate(self.ordered_prof_uids):
            temp_col_dict[prof_uid] = self.prof_pal[i]
        return temp_col_dict

    def _get_ordered_prof_uids(self):
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

    def _make_pcoa_df(self):
        df = pd.read_csv(filepath_or_buffer=self.pcoa_path, sep=',', header=0, index_col=0, skipfooter=1)
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

        return df

    def plot_ordination_figure(self):
        self.large_map_ax.set_extent(extents=(33.0, 41.0, 22.0, 30.0))
        land_110m, ocean_110m, boundary_110m = self._get_naural_earth_features_big_map()
        self._draw_natural_earth_features_big_map(land_110m, ocean_110m,boundary_110m)
        self._put_gridlines_on_large_map_ax()
        self._annotate_big_map()

        # self.zoom_map_ax.set_extent(extents=(38.75, 39.25, 22, 22.5))
        # land_10m, ocean_10m = self._get_naural_earth_features_zoom_map()
        # self._draw_natural_earth_features_zoom_map(land_10m, ocean_10m)
        # self._put_gridlines_on_zoom_map_ax()
        # self._annotate_zoom_map()

        x_list = []
        y_list = []
        color_list = []
        marker_list = []
        for i, sample_uid in enumerate(self.pcoa_df.index.values.tolist()):
            x_list.append(self.pcoa_df['PC1'][sample_uid])
            y_list.append(self.pcoa_df['PC2'][sample_uid])
            sample_name = self.sample_uid_to_sample_name_dict[sample_uid]
            site = self.meta_df.loc[sample_name, 'site']
            site_color = self.site_color_dict[site]
            color_list.append(site_color)
            marker_list.append(self.site_marker_dict[site])

        for x, y, c, m in zip(x_list, y_list, color_list, marker_list):
            self.pc1_pc2_ax.scatter(x, y, c=c, marker=m, edgecolors='black')
        self.pc1_pc2_ax.set_xticks([])
        self.pc1_pc2_ax.set_yticks([])
        self.pc1_pc2_ax.set_xlabel('PC1')
        self.pc1_pc2_ax.set_ylabel('PC2')
        apples = 'asdf'

        y_list = []

        for i, sample_uid in enumerate(self.pcoa_df.index.values.tolist()):
            y_list.append(self.pcoa_df['PC3'][sample_uid])
        for x, y, c, m in zip(x_list, y_list, color_list, marker_list):
            self.pc1_pc3_ax.scatter(x, y, c=c, marker=m, edgecolors='black')
        self.pc1_pc3_ax.set_xticks([])
        self.pc1_pc3_ax.set_yticks([])
        self.pc1_pc3_ax.set_xlabel('PC1')
        self.pc1_pc3_ax.set_ylabel('PC3')

        apples = 'asdf'


        sample_plotting_order_matrix = [[[] for _ in range(len(self.sites))] for _ in range(len(self.ordered_prof_uids))]
        # We will order the samples by mostabund type profile and then by site
        # ugly but fastest is just to go through the df multiple times
        for sample_uid in self.prof_rel_df.index.values.tolist():
            prof_name_ind = self.ordered_prof_uids.index(self.prof_rel_df.loc[sample_uid].idxmax())
            site_ind = self.sites.index(self.meta_df.at[self.sample_uid_to_sample_name_dict[sample_uid], 'site'])
            sample_plotting_order_matrix[prof_name_ind][site_ind].append(sample_uid)
        sample_order = []
        for i in range(len(sample_plotting_order_matrix)):
            for j in range(len(sample_plotting_order_matrix[i])):
                sample_order.extend(sample_plotting_order_matrix[i][j])

        # now lets start by plotting the first 50 and then the remainders
        num_sampls_first_plot = 55
        width = 1 / num_sampls_first_plot
        self._plot_seqs_on_ax(
            ordered_sample_list=sample_order[:num_sampls_first_plot],
            ax = self.sub_plot_seqs_axarr[0],
            width=width,
            x_ind_list = [i * width for i in range(num_sampls_first_plot)])

        self._plot_seqs_on_ax(
            ordered_sample_list=sample_order[num_sampls_first_plot:],
            ax=self.sub_plot_seqs_axarr[1],
            width=width,
            x_ind_list=[i * width for i in range(len(sample_order)-num_sampls_first_plot)])

        from fastkml import kml
        with open(os.path.join(self.input_base_dir, 'test_path.kml'), 'rb') as f:
            doc=f.read()

        k = kml.KML()
        k.from_string(doc)
        print(k.to_string(prettyprint=True))
        apples = 'asdf'

    def _plot_profs_on_ax(self, ordered_sample_list, ax, width, x_ind_list):
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.03, 1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.plot([x_ind_list[0], x_ind_list[-1]+width], [0,0], 'k-', linewidth=ax.spines['right']._linewidth)

        patches_list = []
        color_list = []
        for sample_uid, x_ind in zip(ordered_sample_list, x_ind_list):
            # for each sample we will start at 0 for the y and then add the height of each bar to this
            bottom = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            current_sample_series = self.prof_rel_df.loc[sample_uid]
            non_zero_indices = current_sample_series.nonzero()[0]
            non_zero_sample_series = current_sample_series.iloc[non_zero_indices]
            sample_total = non_zero_sample_series.sum()
            for ser_index, rel_abund in non_zero_sample_series.iteritems():
                height = rel_abund / sample_total
                patches_list.append(Rectangle(
                    (x_ind, bottom), width,
                    height, color=self.prof_color_dict[ser_index]))
                bottom += height
                color_list.append(self.prof_color_dict[ser_index])
        listed_colour_map = ListedColormap(color_list)
        patches_collection = PatchCollection(patches_list, cmap=listed_colour_map)
        patches_collection.set_array(np.arange(len(patches_list)))
        ax.add_collection(patches_collection)

    def _plot_seqs_on_ax(self, ordered_sample_list, ax, width, x_ind_list):
        ax.set_xlim(0, 1)
        prof_depth = 0.3
        ax.set_ylim(-(prof_depth+0.2), 1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])


        ax.plot([x_ind_list[0], x_ind_list[-1]+width], [0,0], 'k-', linewidth=ax.spines['right']._linewidth)
        ax.text(x=-0.05, y=0.5, s='ITS2 sequences', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize='xx-small')
        ax.text(x=-0.05, y=-(prof_depth+0.2)/2, s='predicted\nprofile', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize='xx-small')

        patches_list = []
        color_list = []
        for sample_uid, x_ind in zip(ordered_sample_list, x_ind_list):
            # for each sample we will start at 0 for the y and then add the height of each bar to this
            site = self.meta_df.loc[self.sample_uid_to_sample_name_dict[sample_uid], 'site']
            if site == 'protected':
                ax.scatter(x_ind + (0.5 * width), -0.45, marker=self.site_marker_dict[site],
                           facecolor='black')
            else:
                ax.scatter(x_ind + (0.5*width), -0.45, marker=self.site_marker_dict[site], edgecolor='black', facecolor='white')
            bottom = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            # make the sequence rectangles
            current_sample_series = self.seq_rel_df.loc[sample_uid]
            non_zero_indices = current_sample_series.nonzero()[0]
            non_zero_sample_series = current_sample_series.iloc[non_zero_indices]
            sample_total = non_zero_sample_series.sum()
            for ser_index, rel_abund in non_zero_sample_series.iteritems():
                height = rel_abund / sample_total
                patches_list.append(Rectangle(
                    (x_ind, bottom), width,
                    height, color=self.seq_color_dict[ser_index]))
                bottom += height
                color_list.append(self.seq_color_dict[ser_index])

            # make the profile rectangle
            current_sample_series = self.prof_rel_df.loc[sample_uid]
            non_zero_indices = current_sample_series.nonzero()[0]
            non_zero_sample_series = current_sample_series.iloc[non_zero_indices]
            sample_total = non_zero_sample_series.sum()
            bottom = 0
            for ser_index, rel_abund in non_zero_sample_series.iteritems():
                height = (rel_abund / sample_total)*-prof_depth
                patches_list.append(Rectangle(
                    (x_ind, bottom), width,
                    height, color=self.prof_color_dict[ser_index]))
                bottom += height
                color_list.append(self.prof_color_dict[ser_index])


        listed_colour_map = ListedColormap(color_list)
        patches_collection = PatchCollection(patches_list, cmap=listed_colour_map)
        patches_collection.set_array(np.arange(len(patches_list)))
        ax.add_collection(patches_collection)
        # ax.autoscale_view()
        # ax.figure.canvas.draw()


    def _get_naural_earth_features_big_map(self):
        land_110m = cartopy.feature.NaturalEarthFeature(category='physical', name='land',
                                                       scale='110m')
        ocean_110m = cartopy.feature.NaturalEarthFeature(category='physical', name='ocean',
                                                        scale='110m')
        boundary_110m = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                            name='admin_0_boundary_lines_land', scale='110m')
        return land_110m, ocean_110m, boundary_110m

    def _get_naural_earth_features_zoom_map(self):
        land_10m = cartopy.feature.NaturalEarthFeature(category='physical', name='land',
                                                       scale='50m')
        ocean_10m = cartopy.feature.NaturalEarthFeature(category='physical', name='ocean',
                                                        scale='50m')

        return land_10m, ocean_10m

    def _draw_natural_earth_features_big_map(self, land_110m, ocean_110m, boundary_110m):
        """NB the RGB must be a tuple in a list and the R, G, B must be given as a value between 0 and 1"""
        self.large_map_ax.add_feature(land_110m, facecolor=[(238 / 255, 239 / 255, 219 / 255)],
                                      edgecolor='black')
        self.large_map_ax.add_feature(ocean_110m, facecolor=[(136 / 255, 182 / 255, 224 / 255)],
                                      edgecolor='black')
        self.large_map_ax.add_feature(boundary_110m, edgecolor='gray', linewidth=0.2, facecolor='None')

    def _draw_natural_earth_features_zoom_map(self, land_10m, ocean_10m):
        """NB the RGB must be a tuple in a list and the R, G, B must be given as a value between 0 and 1"""
        self.zoom_map_ax.add_feature(land_10m, facecolor=[(238 / 255, 239 / 255, 219 / 255)],
                                      edgecolor='black')
        self.zoom_map_ax.add_feature(ocean_10m, facecolor=[(136 / 255, 182 / 255, 224 / 255)],
                                      edgecolor='black')


    def _put_gridlines_on_large_map_ax(self):
        """ Although there is a GeoAxis.gridlines() method, this method does not yet allow a lot of
        bespoke options. If we want to only put the labels on the top and left then we have to
        generate a Gridliner object (normally returned by GeoAxis.gridlines() ourselves. We then need
        to manually change the xlabels_bottom and ylabels_right attributes of this Gridliner object.
        We then draw it by adding it to the GeoAxis._gridliners list."""

        xlocs = [32.0, 34.0, 36.0, 38.0, 40.0, 42.0]
        ylocs = [21.0, 23.0, 25.0, 27.0, 29.0, 31.0]

        if xlocs is not None and not isinstance(xlocs, mticker.Locator):
            xlocs = mticker.FixedLocator(xlocs)
        if ylocs is not None and not isinstance(ylocs, mticker.Locator):
            ylocs = mticker.FixedLocator(ylocs)
        g1 = Gridliner(
            axes=self.large_map_ax, crs=ccrs.PlateCarree(), draw_labels=True,
            xlocator=xlocs, ylocator=ylocs)
        g1.xlabels_bottom = False
        g1.ylabels_right = False
        self.large_map_ax._gridliners.append(g1)
        # self.large_map_ax.gridlines(draw_labels=True)

    def _put_gridlines_on_zoom_map_ax(self):
        """ Although there is a GeoAxis.gridlines() method, this method does not yet allow a lot of
        bespoke options. If we want to only put the labels on the top and left then we have to
        generate a Gridliner object (normally returned by GeoAxis.gridlines() ourselves. We then need
        to manually change the xlabels_bottom and ylabels_right attributes of this Gridliner object.
        We then draw it by adding it to the GeoAxis._gridliners list."""

        xlocs = [38.75, 39.0, 39.25, 39.5]
        ylocs = [22.0, 22.25, 22.5, 22.75]

        if xlocs is not None and not isinstance(xlocs, mticker.Locator):
            xlocs = mticker.FixedLocator(xlocs)
        if ylocs is not None and not isinstance(ylocs, mticker.Locator):
            ylocs = mticker.FixedLocator(ylocs)
        g1 = Gridliner(
            axes=self.zoom_map_ax, crs=ccrs.PlateCarree(), draw_labels=True,
            xlocator=xlocs, ylocator=ylocs)
        g1.xlabels_bottom = False
        g1.ylabels_right = False
        self.zoom_map_ax._gridliners.append(g1)
        # self.zoom_map_ax.gridlines(draw_labels=True)

    def _annotate_big_map(self):
        # collect unique tuples of locations
        x_site_coords = [34.934402, 38.960234,39.05165, 39.04878]
        y_site_coords =  [29.514673, 22.303411, 22.26302, 22.26189]

        # coord_sets = set()
        # for i, sample_name in enumerate(self.meta_df.index.values.tolist()):
        #     coord_sets.add((self.meta_df['collection_longitude'][sample_name], self.meta_df['collection_latitude'][sample_name]))
        # x_site_coords = []
        # y_site_coords = []
        # for site in coord_sets:
        #     x_site_coords.append(site[0])
        #     y_site_coords.append(site[1])

        self.large_map_ax.plot(x_site_coords[0], y_site_coords[0], 'o', markerfacecolor='white', markeredgecolor='black', markersize=8)
        self.large_map_ax.plot(x_site_coords[1], y_site_coords[1], '^', markerfacecolor='black', markeredgecolor='black', markersize=8)
        self.large_map_ax.plot(x_site_coords[2], y_site_coords[2], '^', markerfacecolor='gray', markeredgecolor='black', markersize=8)
        self.large_map_ax.plot(x_site_coords[3], y_site_coords[3], '^', markerfacecolor='white', markeredgecolor='black', markersize=8)

    def _annotate_zoom_map(self):
        # collect unique tuples of locations
        x_site_coords = [ 39.05165, 39.04878, 38.960234]
        y_site_coords = [ 22.26302, 22.26189, 22.303411]

        self.zoom_map_ax.plot(x_site_coords[0], y_site_coords[0], 'bo')
        self.zoom_map_ax.plot(x_site_coords[1], y_site_coords[1], 'go')
        self.zoom_map_ax.plot(x_site_coords[2], y_site_coords[2], 'ko')

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
sof = SampleOrdinationFigure()
sof.plot_ordination_figure()