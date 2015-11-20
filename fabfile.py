from fabric.api import run, cd, hosts, env

env.use_ssh_config = True
env.hosts = {'galaxy': ['galaxy']}

@hosts('galaxy')
def start():
    with cd('/home/galaxy/galaxy-dist'):
        run('GALAXY_RUN_ALL=1 sh run.sh --daemon')

@hosts('galaxy')
def restart():
    with cd('/home/galaxy/galaxy-dist'):
        run('GALAXY_RUN_ALL=1 sh run.sh --stop-daemon')
        run('GALAXY_RUN_ALL=1 sh run.sh --daemon')

@hosts('galaxy')
def deploy_vm():
    with cd('/home/galaxy/geec_galaxy'):
        run('git pull')

