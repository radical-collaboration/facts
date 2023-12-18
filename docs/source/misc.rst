.. _chapter_limitations:

Known Limitations
=========

``tlm/sterodynamics`` relies upon the CMIP6 archive to learn the correlation between global mean
thermosteric sea-level rise and ocean dynamic sea level. It does not currently have a way to use
a separate label for the scenario used for global ocean heat content (which uses output from
the ``fair/temperature`` module) and the scenario used for this correlation. Accordingly, it will
not currently work if ``scenario`` is specified to be anything other than either (1) one of the
five ScenarioMIP SSPs used by AR6 (ssp119, ssp126, ssp245, ssp370, ssp585) or (2) a warming
level-based scenario with a name given by tlimX.XwinY.Y where X.X specified the warming level
target in 2081-2100 and Y.Y the width of the target. If you are trying to run a scenario that
does not fit one of these patterns, you should rename it so that it does. A more general capability will
be added in a future release.

``extremesealevel/pointsoverthreshold`` relies upon GESLA tide-gauge data for estimating
the historical extreme sea level distribution. It searches for the geographically nearest
tide gauge corresponding to a particular local sea level projection. Thus, runs producing extreme sea level output should be limited to
an intentionally chosen, targeted set of tide gauge locations.

The standard pipeline for ``extremesealevel/pointsoverthreshold``
uses a trimmed data file, extremesealevel_pointsoverthreshold_data.tgz.
This is because rhe full GESLA data set (extremesealevel_pointsoverthreshold_fulldata.tgz)
is quite large, and it is particularly slow because
it is broken into many files. You may either want to adjust
the pipeline to use the full data file or produce
your own trimmed data file that includes the GESLA sites
of interest to you. 

On the Rutgers Amarel system, some modules --
especially emulandice -- exhibit hard-to-trace problems that appear to 
be sensitive to the version of Python used. On Amarel, runs have been 
successfully completed with Python 3.9.6 and 3.10.12, but not 3.6.8. In 
the Docker container, 3.8.10 is know to run successfully. If you are
having hard-to-identify run issues, consider checking the version of Python you are using. 

