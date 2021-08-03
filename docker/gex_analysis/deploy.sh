
# Prereq: build.sh and test.sh run
REGISTRY=$1
ORG=$2
VERSION=$3
USER=$3
PSSWD=$4

echo "\n LOGIN \n"
docker login $REGISTRY -u $USER --password $PSSWD

echo "\n PREPARING \n"

echo "\n PUSHING \n"

docker push $REGISTRY/$ORG/gex_analysis:$VERSION