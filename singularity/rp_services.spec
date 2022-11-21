
Bootstrap: library
From: ubuntu:22.04

%environment
    export LC_ALL=C

%post
    chmod 1777 /tmp
    mkdir -p   /tmp/mongodb /tmp/rabbitmq
    chmod 0755 /tmp/mongodb /tmp/rabbitmq
    . /tmp/rp_services.env
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
    . /tmp/rp_services.env
    mongod --config /tmp/mongodb.conf
    env | grep RABBITMQ
    rabbitmq-server


