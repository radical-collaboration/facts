import pandas as pd
import xarray as xr
import warnings

#python script downloaded from https://gesla787883612.wordpress.com/downloads/
#TH: adapted to fix some bugs

class GeslaDataset:
    """A class for loading data from GESLA text files into convenient in-memory
    data objects. By default, single file requests are loaded into
    `pandas.DataFrame` objects, which are similar to in-memory spreadsheets.
    Multifile requests are loaded into `xarray.Dataset` objects, which are
    similar to in-memory NetCDF files."""

    def __init__(self, meta_file, data_path):
        """Initialize loading data from a GESLA database.

        Args:
            meta_file (string): path to the metadata file in .csv format.
            data_path (string): path to the directory containing GESLA data
                files.
        """
        self.meta = pd.read_csv(meta_file)
        self.meta.columns = [
            c.replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("/", "_")
            .lower()
            for c in self.meta.columns
        ]
        self.meta.loc[:, "start_date_time"] = [
            pd.to_datetime(d,format='mixed') for d in self.meta.loc[:, "start_date_time"]
        ]
        self.meta.loc[:, "end_date_time"] = [
            pd.to_datetime(d,format='mixed') for d in self.meta.loc[:, "end_date_time"]
        ]
        self.data_path = data_path
        self.meta["filename"] = self.construct_filenames()

    def construct_filenames(self):
        return self.meta.apply(
            lambda x: x.loc["site_name"].lower()
            + "-"
            + x.loc["site_code"].lower()
            + "-"
            + x.loc["country"].lower()
            + "-"
            + x.loc["contributor_abbreviated"].lower(),
            axis=1,
        )

    def file_to_pandas(self, filename, return_meta=True):
        """Read a GESLA data file into a pandas.DataFrame object. Metadata is
        returned as a pandas.Series object.

        Args:
            filename (string): name of the GESLA data file. Do not prepend path.
            return_meta (bool, optional): determines if metadata is returned as
                a second function output. Defaults to True.

        Returns:
            pandas.DataFrame: sea-level values and flags with datetime index.
            pandas.Series: record metadata. This return can be excluded by
                setting return_meta=False.
        """
        with open(self.data_path + filename, "r") as f:
            data = pd.read_csv(
                f,
                skiprows=41,
                names=["date", "time", "sea_level", "qc_flag", "use_flag"],
                sep="\s+",
                parse_dates=[[0, 1]],
                index_col=0,
            )
            if data.index[data.index.duplicated()].size > 0:
                #data = data.drop_duplicates() #this removes any duplicate row, rather than duplicate timestamps
                data = data[~data.index.duplicated(keep='first')]
                
            if return_meta:
                try:
                    meta = self.meta.loc[self.meta['file_name'] == filename].iloc[0]
                except:
                    meta = self.meta.loc[self.meta.filename == filename].iloc[0]
                #print(self.meta.filename[2495])
                #meta = self.meta.loc[self.meta.filename == filename].iloc[0]
                return data, meta
            else:
                return data

    def files_to_xarray(self, filenames):
        """Read a list of GESLA filenames into a xarray.Dataset object. The
        dataset includes variables containing metadata for each record.

        Args:
            filenames (list): list of filename strings.

        Returns:
            xarray.Dataset: data, flags, and metadata for each record.
        """
        data = xr.concat(
            [
                self.file_to_pandas(f, return_meta=False).to_xarray()
                for f in filenames
            ],
            dim="station",
        )

        idx = [
            s.Index for s in self.meta.itertuples() if s.filename in filenames
        ]
        meta = self.meta.loc[idx]
        meta.index = range(meta.index.size)
        meta.index.name = "station"
        data = data.assign({c: meta[c] for c in meta.columns})

        return data

    def load_N_closest(self, lat, lon, N=1, force_xarray=False):
        """Load the N closest GESLA records to a lat/lon location into a
        xarray.Dataset object. The dataset includes variables containing
        metadata for each record.

        Args:
            lat (float): latitude on the interval [-90, 90]
            lon (float): longitude on the interval [-180, 180]
            N (int, optional): number of locations to load. Defaults to 1.
            force_xarray (bool, optional): if N=1, the default behavior is to
                return a pandas.DataFrame object containing data/flags and a
                pandas.Series object containing metadata. Set this argument to
                True to return a xarray Dataset even if N=1. Defaults to False.

        Returns:
            xarray.Dataset: data, flags, and metadata for each record.
        """
        N = int(N)
        if N <= 0:
            raise Exception("Must specify N > 0")

        d = (self.meta.longitude - lon) ** 2 + (self.meta.latitude - lat) ** 2
        idx = d.sort_values().iloc[:N].index
        meta = self.meta.loc[idx]

        if (N > 1) or force_xarray:
            return self.files_to_xarray(meta.filename.tolist())

        else:
            data, meta = self.file_to_pandas(meta.filename.values[0])
            return data, meta

    def load_lat_lon_range(
        self,
        south_lat=-90,
        north_lat=90,
        west_lon=-180,
        east_lon=180,
        force_xarray=False,
    ):
        """Load GESLA records within a rectangular lat/lon range into a xarray.
        Dataset object.

        Args:
            south_lat (float, optional): southern extent of the range. Defaults
                to -90.
            north_lat (float, optional): northern extent of the range. Defaults
                to 90.
            west_lon (float, optional): western extent of the range. Defaults
                to -180.
            east_lon (float, optional): eastern extent of the range. Defaults
                to 180.
            force_xarray (bool, optional): if there is only one record in the
                lat/lon range, the default behavior is to return a
                pandas.DataFrame object containing data/flags and a
                pandas.Series object containing metadata. Set this argument to
                True to return a xarray.Dataset even if only one record is
                selected. Defaults to False.

        Returns:
            xarray.Dataset: data, flags, and metadata for each record.
        """
        if west_lon > 0 & east_lon < 0:
            lon_bool = (self.meta.longitude >= west_lon) | (
                self.meta.longitude <= east_lon
            )
        else:
            lon_bool = (self.meta.longitude >= west_lon) & (
                self.meta.longitude <= east_lon
            )
        lat_bool = (self.meta.latitude >= south_lat) & (
            self.meta.latitude <= north_lat
        )
        meta = self.meta.loc[lon_bool & lat_bool]

        if (meta.index.size > 1) or force_xarray:
            return self.files_to_xarray(meta.filename.tolist())

        else:
            data, meta = self.file_to_pandas(meta.filename.values[0])
            return data, meta