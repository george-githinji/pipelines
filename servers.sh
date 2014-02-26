tmux new-session -s servers -d "ssh ssh.sanger.ac.uk"
tmux split-window -v "ssh ggithinji@172.16.12.144"
# tmux split-window -v "ssh labserver"
tmux attach -t servers
