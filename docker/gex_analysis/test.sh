
# Prereq: build.sh run

REGISTRY=$1
ORG=$2

echo "\n PREPARING TO TEST \n"

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
echo "Getting tag $TAG"

docker exec -t $REGISTRY/$ORG/gex_analysis:$TAG pytest --pyargs gex_analysis
