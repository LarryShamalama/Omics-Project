Here is the command to create an RSA key
`
ssh-keygen -t rsa -f ./key -C larry_dong_10
# password for key is required
chmod 400 key
cat key.pub # see contents of key.pub
`
And copy paste the contents into the Google Cloud virtual environment

:exclamation: Posting a key for the public is a dangerous and bad habit
