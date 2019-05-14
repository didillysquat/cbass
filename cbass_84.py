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
        self.profile_uid_to_profile_name_dict = {}
        self.profile_name_to_profile_uid_dict = {}
        self.prof_rel_df = self._make_prof_rel_df()
        self.fig = plt.figure(figsize=(8, 12))
        self.gs = gridspec.GridSpec(2, 2)
        self.large_map_ax = plt.subplot(self.gs[:1, :1], projection=ccrs.PlateCarree(), zorder=1)
        # self.zoom_map_ax = plt.subplot(self.gs[:1, 1:2], projection=ccrs.PlateCarree(), zorder=1)
        self.pc1_pc2_ax = plt.subplot(self.gs[1:2, :1])
        self.site_color_dict = {'eilat':'white', 'kaust_0': 'black', 'kaust_1': 'gray', 'kaust_2':'white'}
        self.pc1_pc3_ax = plt.subplot(self.gs[1:2, 1:2])
        self.sub_plot_gs = gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec=self.gs[0:1, 1:2])
        self.sub_plot_seqs = plt.subplot(self.sub_plot_gs[0:1,0:1])
        self.sub_plot_profiles = plt.subplot(self.sub_plot_gs[1:2, 0:1])
        self.site_marker_dict = {'eilat': 'o', 'kaust_0': '^', 'kaust_1': '^', 'kaust_2': '^' }

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
        return df

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

        return df

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
                site_name.append('kaust_0')
            elif df['collection_latitude'][sample_name] == 22.26189:
                site_name.append('kaust_1')
            elif df['collection_latitude'][sample_name] == 22.303411:
                site_name.append('kaust_2')
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



        self.sub_plot_seqs = plt.subplot(self.sub_plot_gs[0:1, 0:1])
        self.sub_plot_profiles = plt.subplot(self.sub_plot_gs[0:1, 0:1])




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
        x_site_coords = [34.934402, 39.05165, 39.04878, 38.960234]
        y_site_coords =  [29.514673, 22.26302, 22.26189, 22.303411]

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

sof = SampleOrdinationFigure()
sof.plot_ordination_figure()