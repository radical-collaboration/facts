#!/bin/sh

tmp="$(pwd)/rct_tmp/"
uid="$(id -un)"

singularity instance list \
    | grep -v INSTANCE \
    | grep rct_ \
    | cut -f 1 -d ' ' \
    | xargs -r -n 1 singularity instance stop

test -s "$tmp/rct_mongodb.pid"  && kill -9 $(cat "$tmp/rct_mongodb.pid")
test -s "$tmp/rct_rabbitmq.pid" && kill -9 $(cat "$tmp/rct_rabbitmq.pid")

ps -f -u $uid \
    | grep -v grep \
    | grep -e '/erlang/' \
    | awk '{print $2}' \
    | xargs -r kill -9

ps -f -u $uid \
    | grep -v grep \
    | grep -e '\<mongod\>' \
    | awk '{print $2}' \
    | xargs -r kill -9

rm -f "$tmp/rct_mongodb.pid"
rm -f "$tmp/rct_rabbitmq.pid"

