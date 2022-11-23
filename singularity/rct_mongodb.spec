
Bootstrap: library
From: ubuntu:22.04

%environment
    export LC_ALL=C

%post
    chmod 1777 /tmp
    mkdir -p   /tmp/mongodb
    chmod 0755 /tmp/mongodb
    apt update
    apt install -y gpg wget
    wget -qO - https://www.mongodb.org/static/pgp/server-6.0.asc \
        | gpg --dearmor \
        >> /usr/share/keyrings/mongodb.gpg
    echo "deb [ arch=amd64,arm64 signed-by=/usr/share/keyrings/mongodb.gpg ]" \
         "https://repo.mongodb.org/apt/ubuntu jammy/mongodb-org/6.0 multiverse"\
        >> /etc/apt/sources.list.d/mongodb-org-6.0.list
    apt update
    apt install -y mongodb-org rabbitmq-server

%runscript
    mongod --config /tmp/mongodb.conf &
    echo "$!" > /tmp/mongodb.pid

