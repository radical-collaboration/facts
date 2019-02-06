# Instructions to run EnTK 0.7 scripts

## Installation and environment setup on laptop/VM

* Virtualenv creation

```
virtualenv $HOME/myenv
source $HOME/myenv/bin/activate
```

* EnTK Installation: You can directly install from pypi

```
pip install radical.entk
```

You can check the version with ```radical-stack```, this
should print 0.7.* (where * is a number indicating versions
of patches).

NOTE: If your target machine is Blue Waters, note the radical
pilot version printed.

* RabbitMQ

A docker instance is already setup on radical.two with
hostname 'two.radical-project.org' and port '33239' that
is specified in the resource config file

If you want to create your own, please see [EnTK docs](https://radicalentk.readthedocs.io/en/latest/install.html#installing-rabbitmq).

* MongoDB

You can create a new mongodb instance on [mlab](https://mlab.com/).
Steps broadly include:

1. Create a new db instance with approriate geographic location,
instance size, etc.
2. Once created, select the db and add an user to it. Specify
a user name and password.
3. Copy the db url and specify it as environment variable

```
export RADICAL_PILOT_DBURL='<mongo url>'
```

## Executing your script

To execute your script, please first set the verbosity level:
```
export RADICAL_ENTK_VERBOSE=REPORT
```
We recommend setting it INFO first. Although it generates a lot more
messages, the first time you are running your script probably will
require debugging at some point.

You can execute your script with:
```
python runme.py
```

## Issue reporting

If there are any issues/questions, please create a ticket in this 
Github repository. After investigating the issue, I can move it to
EnTK/RP/SAGA.
