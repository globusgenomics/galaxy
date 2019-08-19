## requires jq
# get files by id
output=$(curl -L -H "x-access-token: $1" https://api.basespace.illumina.com/v1pre3/files/$2 | jq '.Response.Name' -r)
echo $output
curl -L -H "x-access-token: $1" https://api.basespace.illumina.com/v1pre3/files/$2/content > $3
