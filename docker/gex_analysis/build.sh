REGISTRY=$1
ORG=$2

echo "\n PREPARING \n"

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
echo "Getting tag $TAG"

echo "\n BUILDING \n"

docker build -t $REGISTRY/$ORG/gex_analysis:$TAG . --no-cache

docker push $REGISTRY/$ORG/gex_analysis:$TAG
