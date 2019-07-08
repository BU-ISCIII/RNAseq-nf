Bootstrap: docker
From: buisciii/centos7_base_image:latest

%files
    ./scif_app_recipes/ /opt/
%post
    echo "Install basic development tools"
    yum -y groupinstall "Development Tools"
    yum -y update && yum -y install wget curl openssl-devel geos-devel udunits2-devel libxml2-devel cairo-devel libgit2-devel


    echo "Install python2.7 setuptools and pip"
    yum -y install python-setuptools
    easy_install pip


    echo "Install miniconda"
    curl -fsSL https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -o miniconda_v4.6.14.sh
    bash miniconda_v4.6.14.sh -b -p /opt/miniconda
    echo 'export PATH="/opt/miniconda/bin:$PATH"' >> $SINGULARITY_ENVIRONMENT
    ln -s /opt/miniconda/bin/conda bin
    conda init; conda config --append channels conda-forge; conda config --add channels defaults; conda config --add channels bioconda


    echo "Installing SCI-F"
    pip install scif

    echo "Installing FastQC app" && \
    scif install /opt/scif_app_recipes/fastqc_v0.11.7_centos7.scif && \
    echo "Installing trimmomatic app" && \
    scif install /opt/scif_app_recipes/trimmomatic_v0.38_centos7.scif && \
    echo "Installing STAR app" && \
    scif install /opt/scif_app_recipes/STAR_v2.6.1d_centos7.scif       
    echo "Installing Hisat2 app" && \
    scif install /opt/scif_app_recipes/hisat2_v2.1.0_centos7.scif
    echo "Installing Picard app" && \
    scif install /opt/scif_app_recipes/picard_v2.18.27_centos7.scif    
    echo "Installing csvtk app" && \
    scif install /opt/scif_app_recipes/csvtk_v0.17.0_centos7.scif
    echo "Installing preseq app" && \
    scif install /opt/scif_app_recipes/preseq_v2.0.3_centos7.scif
    echo "Installing R app" && \
    scif install /opt/scif_app_recipes/R_v3.5.1_centos7.scif && \
    echo "Installing RSeQC app" && \
    scif install /opt/scif_app_recipes/RSeQC_v3.0.0_centos7.scif
    echo "Installing samtools app" && \
    scif install /opt/scif_app_recipes/samtools_v1.9_centos7.scif && \
    echo "Installing stringtie app" && \
    scif install /opt/scif_app_recipes/stringtie_v1.3.5_centos7.scif && \
    echo "Installing subread app" && \
    scif install /opt/scif_app_recipes/subread_v1.6.4_centos7.scif && \
    echo "Installing gffread app" && \
    scif install /opt/scif_app_recipes/gffread_v0.9.12_centos7.scif && \
    echo "Installing deeptools app" && \
    scif install /opt/scif_app_recipes/deeptools_v2.5.4_centos7.scif && \
    echo "Installing multiqc app" && \
    scif install /opt/scif_app_recipes/multiqc_v1.7_centos7.scif
    
    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It would be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\${PATH}:$in" >> $SINGULARITY_ENVIRONMENT;done

    find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:$in" >> $SINGULARITY_ENVIRONMENT ;done

    if [[ ":$PATH:" == *":/scif/apps/snppipeline:"* ]];then

        export CLASSPATH=/scif/apps/picard/picard.jar:$CLASSPATH >> $SINGULARITY_ENVIRONMENT

    fi
    ## R packages
    
    # Install core R dependencies
    echo "r <- getOption('repos'); r['CRAN'] <- 'https://ftp.acc.umu.se/mirror/CRAN/'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('dupRadar',dependencies=TRUE,lib='/usr/local/lib64/R/library')" && \
    Rscript -e "install.packages('data.table',dependencies=TRUE,lib='/usr/local/lib64/R/library')" && \
    Rscript -e "install.packages('gplots',dependencies=TRUE,lib='/usr/local/lib64/R/library')" && \
    Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('edgeR',dependencies=TRUE,lib='/usr/local/lib64/R/library')" && \
    Rscript -e "install.packages('markdown',dependencies=TRUE,lib='/usr/local/lib64/R/library')"

%runscript
    exec scif "$@"
