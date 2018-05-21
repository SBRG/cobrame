FROM sbrg/cobrame:qminos AS qminos
FROM sbrg/cobrame:dependencies
USER root

# Get qminos shared libraries
COPY --from=qminos /source/libminos.a .
COPY --from=qminos /source/libquadminos.a .

RUN echo "@community http://dl-cdn.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories

# If soplex is present in directory, copy it into source
# nothing.txt is present as a workaround to prevent COPY from returning an
# error. Do not delete it.
ENV SOPLEX_VERSION=3.1.1
# Add qminos solver
COPY nothing.txt soplex-$SOPLEX_VERSION.tg* /source/

# Create user with UID=1000 and in the 'users' group
# and ensure these dirs writable by 'users' group

# Shadow allows use of useradd.
ENV PACKAGES="\
	shadow \
"
RUN apk add --no-cache $PACKAGES

ENV NB_USER=meuser \
	NB_UID=1000 \
	NB_GID=100 \
	HOME=/home/$NB_USER

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
	chmod -R a+rwx /etc/passwd /etc/group /home/ /usr/lib/python$PYTHON_VERSION/site-packages && \
	chmod -R a+rwx /usr/bin /source/

# Switch back to non-root
USER $NB_UID

WORKDIR /home/$NB_USER

RUN echo \
	&& cd /source \
	# Install cobrapy version 0.5.11
	&& git clone https://github.com/opencobra/cobrapy.git \
	&& cd /source/cobrapy \
	&& git checkout 0.5.11 \
	&& $PYTHON setup.py install \
	&& cd /source \
	# if soplex was copied into source, install soplex_cython
	&& if [[ -e /source/soplex-$SOPLEX_VERSION.tgz ]]; then \
        git clone https://github.com/SBRG/soplex_cython.git; \
        cd /source/soplex_cython; \
        mv /source/soplex-$SOPLEX_VERSION.tgz  /source/soplex_cython/ ;\
        pip install . ;\
        fi \

    # Install remaining ME-model software. qMINOS/solvemepy is automatically
    # installed
    && pip install escher \
    && cd /source \
	&& git clone https://github.com/SBRG/cobrame.git \
	&& git clone https://github.com/SBRG/ecolime.git \
	&& git clone https://github.com/SBRG/solvemepy.git \
	&& cd /source/cobrame \
	&& $PYTHON  setup.py develop \
	&& cd /source/solvemepy \
	&& git checkout tags/v1.0.1 \
	&& cp /source/libminos.a ./ \
	&& cp /source/libquadminos.a ./ \
	&& $PYTHON  setup.py develop \
	&& cd /source/ecolime \
	&& $PYTHON  setup.py develop \

	# build iJL1678b ME-model
	&& $PYTHON  /source/ecolime/ecolime/build_me_model.py \
	&& cp -r /source/ecolime/ecolime/me_models /home/$NB_USER/ \
	&& cp /source/ecolime/ecolime/solve_demo.ipynb /home/$NB_USER/ \
	&& echo
