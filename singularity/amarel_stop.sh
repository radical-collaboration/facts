#!/bin/sh

uid="$(id -un)"

singularity instance list \
    | grep -v INSTANCE \
    | grep rct_ \
    | cut -f 1 -d ' ' \
    | xargs -r -n 1 singularity instance stop

ps -f -u $uid \
    | grep -v grep \
    | grep -e '\<mongod\>' \
    | awk '{print $2}' \
    | xargs -r kill -9

ps -f -u $uid \
    | grep -v grep \
    | grep -e '\<erlang\>' \
    | awk '{print $2}' \
    | xargs -r kill -9

ps -f -u $uid \
    | grep -v grep \
    | grep -e '\<beam.smp\>' \
    | awk '{print $2}' \
    | xargs -r kill -9

