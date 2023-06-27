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
