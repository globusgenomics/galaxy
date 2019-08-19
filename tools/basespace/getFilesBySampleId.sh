## requires jq

# get files by sample id
mkdir $3
for fileId in $(curl -L -H "x-access-token: $1" https://api.basespace.illumina.com/v1pre3/samples/$2/files | jq '.Response.Items[].Id' -r)
do
        output=$(curl -L -H "x-access-token: $1" https://api.basespace.illumina.com/v1pre3/files/$fileId | jq '.Response.Name' -r)
        curl -L -H "x-access-token: $1" https://api.basespace.illumina.com/v1pre3/files/$fileId/content > $3/$output;
done

