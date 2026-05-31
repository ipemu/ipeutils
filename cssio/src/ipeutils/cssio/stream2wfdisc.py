
from obspy import read
import datetime
import os
import numpy as np

# CSS3.0 WFDISC format specification:
# https://anf.ucsd.edu/publications/css30.pdf

# initialize the wfdisc data with null values
wfdisc_NULL = {
    'sta': '-',
    'chan': '-',
    'time': -9999999999.999,
    'wfid': -1,
    'chanid': -1,
    'jdate': -1,
    'endtime': 9999999999.999,
    'nsamp': -1,
    'samprate': -1.0,
    'calib': 0.0,
    'calper': -1.0,
    'instype': '-',
    'segtype': '-',
    'datatype': '-',
    'clip': '-',
    'dir': '-',
    'dfile': '-',
    'foff': 0,
    'commid': -1,
    'lddate': '-'
}

def write_wfdisc(wfd, wfdisc_data):
    """
    Write one line of WFDISC data to the specified file.

    :param wfdisc_data: A dictionary containing the WFDISC data.
    :type wfdisc_data: dict
    :param wfd: A file object opened for writing the WFDISC data.
    :type wfd: file object
    """
    # Define the format string for writing the data (283 characters total)
    format_string = "%-6.6s %-8.8s %17.5f %8ld %8ld %8ld %17.5f %8ld %11.7f %16.6f %16.6f %-6.6s %1s %-2.2s %1s %-64.64s %-32.32s %10ld %8ld %-17.17s\n"
    
    # Write the data to the file using the format string
    #with open(filename, 'at') as f:
    wfd.write(format_string % (
            wfdisc_data['sta'],
            wfdisc_data['chan'],
            wfdisc_data['time'],
            wfdisc_data['wfid'],
            wfdisc_data['chanid'],
            wfdisc_data['jdate'],
            wfdisc_data['endtime'],
            wfdisc_data['nsamp'],
            wfdisc_data['samprate'],
            wfdisc_data['calib'],
            wfdisc_data['calper'],
            wfdisc_data['instype'],
            wfdisc_data['segtype'],
            wfdisc_data['datatype'],
            wfdisc_data['clip'],
            wfdisc_data['dir'],
            wfdisc_data['dfile'],
            wfdisc_data['foff'],
            wfdisc_data['commid'],
            wfdisc_data['lddate']
        ))

def tr_stats_to_wfdisc(tr_stats,wf,wfid=1,chanid=-1):
    """
    Convert an obspy Trace.stats to CSS3.0 WFDISC format.

    :type tr_stats: obspy.Trace.stats
    :param tr_stats: The stats object from an obspy Trace.
    :type wf: binary file object
    :param wf: file for writing the waveform data
    :param wfid: The incremental waveform id for this Trace.
    :type wfid: int
    :param chanid: The channel id for this Trace
            (default -1 indicates unknown channel id).
    :type chanid: int
    :rtype: dict
    :return: A dictionary containing the WFDISC data.
    """
    wfdisc_data = wfdisc_NULL.copy()
    wfdisc_data['sta'] = tr_stats.station
    wfdisc_data['chan'] = tr_stats.channel
    wfdisc_data['time'] = tr_stats.starttime.timestamp
    wfdisc_data['wfid'] = wfid # incremental waveform id
    wfdisc_data['chanid'] = chanid # channel id
    wfdisc_data['jdate'] = int(tr_stats.starttime.julday)
    wfdisc_data['endtime'] = tr_stats.endtime.timestamp
    wfdisc_data['nsamp'] = tr_stats.npts
    wfdisc_data['samprate'] = tr_stats.sampling_rate
    wfdisc_data['calib'] = tr_stats.calib if hasattr(tr_stats, 'calib') else 1.0
    wfdisc_data['calper'] = tr_stats.calper if hasattr(tr_stats, 'calper') else -1.0
    wfdisc_data['instype'] = '-' # unknown instrument type
    wfdisc_data['segtype'] = 'V' # velocity data segment type 
    wfdisc_data['datatype'] = 'f4' # little-endian float32 format for waveform data 
    wfdisc_data['clip'] = '-' # default clipped flag
    wfdisc_data['dir'] = '.' # default directory
    wfdisc_data['dfile'] = os.path.basename(wf.name) 
    wfdisc_data['foff'] = wf.tell() # the current file offset for the waveform data 
    wfdisc_data['commid'] = -1 # default comment id
    # actual now UTC date and time in format YYYY-MM-DD HH:MM:SS
    wfdisc_data['lddate'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M')
    return wfdisc_data

def save_stream_to_wfdisc(st,filename_prefix='output'):
    """
    Save an obspy Stream to one wfdisc text file and one waveform binary file in CSS3.0 format.

    :param st: The obspy Stream to save.
    :type st: obspy.Stream
    :param filename_prefix: The prefix for the output files. Including path if needed.
    :type filename_prefix: str

    Note: The function will create two files:
    1. A text file with the extension ".wfdisc" containing the WFDISC metadata
       for each Trace in the Stream.
    2. A binary file with the extension ".w" containing the waveform data
       for all Traces in the Stream in CSS3.0 format (f4).
       Files are opened in write mode, so existing files with the same names
       will be overwritten.
       The WFDISC metadata is written in a fixed-width format
       according to the CSS3.0 specification,
       and the waveform data is written in little-endian float32 format.
    """
    # remove ".wfdisc" extension if it exists in the filename prefix
    if filename_prefix.endswith('.wfdisc'):
        filename_prefix = filename_prefix[:-8]
    # create names of the output files
    file = filename_prefix + '.wfdisc'
    file_bin = filename_prefix + '.w'
    wfd = open(file, 'w') # open the wfdisc text file for writing
    # open the binary file for writing the waveform data
    wf = open(file_bin, 'wb')
    wfid = 0
    for tr in st:
        wfid += 1
        wfdisc_data = tr_stats_to_wfdisc(tr.stats,wf,wfid)
        write_wfdisc(wfd, wfdisc_data)
        # write the waveform data to the binary file in CSS3.0 format (f4)
        # convert the data to float32 format
        data = tr.data.astype(np.float32)
        # write the data to the binary file
        data.tofile(wf)
    wf.close()

# Example usage:
if __name__ == "__main__":
    # Read a sample seismic data file using obspy
    st = read()
    # Save the stream to wfdisc format
    save_stream_to_wfdisc(st)
