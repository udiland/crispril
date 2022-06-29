FROM tiangolo/uwsgi-nginx-flask:python3.7

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.490-1_amd64.deb && dpkg -i mafft_7.490-1_amd64.deb

RUN echo 'export PATH="/usr/bin" >> ~/.bashrc'

RUN wget http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz \
	&& tar -zxvf phylip-3.697.tar.gz \
	&& cd ./phylip-3.697/src \
	&& apt update \
	&& apt install gcc -y \
	&& make -f Makefile.unx install 


#RUN echo -e '#!/bin/bash\n echo phylip-3.697/exe/./protdist "$@"' > /usr/bin/protdist && \
#	chmod +x /usr/bin/protdist

RUN mkdir /crispys_out
RUN mkdir /input

ENV PATH="/usr/bin:/crispys_out:${PATH}" 

COPY ./app /app
COPY ./crispys_code /app/crispys_code


ENV PYTHONPATH=/app/crispys_code
RUN pip install -r /app/requirements.txt

RUN cp ./phylip-3.697/exe/protdist /crispys_out/

