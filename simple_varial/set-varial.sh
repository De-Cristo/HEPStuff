python -c "import varial_ext" &> /dev/null
varial_nonexisting=$?

if [ $varial_nonexisting != 0 ]; then
    if [ -f Varial/setup.py ]; then
        echo "Updating Varial."
        cd Varial
#         git pull
        cd -
    else
        echo "Installing Varial."
        git clone git@github.com:HeinerTholen/Varial.git
    fi

    export PYTHONPATH=$PWD/Varial:$PYTHONPATH
    export PATH=$PATH:$PWD/Varial/bin
fi