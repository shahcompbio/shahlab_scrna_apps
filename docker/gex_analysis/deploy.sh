
# Prereq: build.sh and test.sh run
REGISTRY=$1
ORG=$2
USER=$3
PSSWD=$4

echo "\n LOGIN \n"
docker login $REGISTRY -u $USER --password $PSSWD

echo "\n PREPARING \n"

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
echo "Getting tag $TAG"

echo "\n PUSHING \n"

docker push $REGISTRY/$ORG/gex_analysis:$TAG