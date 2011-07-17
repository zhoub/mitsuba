<?xml version='1.0' encoding='utf-8'?>

<!-- Stylesheet to upgrade pre-0.3.0 scenes -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:output method="xml" indent="yes" encoding="utf-8"/>
	<xsl:preserve-space elements="*"/>

	<xsl:template match="scene">
		<scene>
			<xsl:attribute name="version">0.3.0</xsl:attribute>
			<xsl:apply-templates select="node()"/>
		</scene>
	</xsl:template>

	<!-- Replace the old version of the lookAt tag -->
	<xsl:template match="lookAt[@ox]">
		<scale x="-1"/>
		<lookAt>
			<xsl:attribute name="origin">
				<xsl:value-of select="@ox"/>
				<xsl:text>, </xsl:text>
				<xsl:value-of select="@oy"/>
				<xsl:text>, </xsl:text>
				<xsl:value-of select="@oz"/>
			</xsl:attribute>
			<xsl:attribute name="target">
				<xsl:value-of select="@tx"/>
				<xsl:text>, </xsl:text>
				<xsl:value-of select="@ty"/>
				<xsl:text>, </xsl:text>
				<xsl:value-of select="@tz"/>
			</xsl:attribute>
			<xsl:if test="@ux">
				<xsl:attribute name="up">
					<xsl:value-of select="@ux"/>
					<xsl:text>, </xsl:text>
					<xsl:value-of select="@uy"/>
					<xsl:text>, </xsl:text>
					<xsl:value-of select="@uz"/>
				</xsl:attribute>
			</xsl:if>
		</lookAt>
	</xsl:template>

	<!-- There are no more 'diffuseAmount' or 'specularAmount' parameters in
		 the microfacet/phong/ward plugins, and the default values have changed -->
	<xsl:template match="bsdf[@type='microfacet' or @type='phong' or @type='ward']">
		<xsl:variable name="diffuseAmount">
			<xsl:choose>
				<xsl:when test="float[@name='diffuseAmount']"> 
					<xsl:value-of select="float[@name='diffuseAmount']/@value"/>
				</xsl:when>
				<xsl:otherwise>1.0</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="specularAmount">
			<xsl:choose>
				<xsl:when test="float[@name='specularAmount']"> 
					<xsl:value-of select="float[@name='specularAmount']/@value"/>
				</xsl:when>
				<xsl:otherwise>1.0</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="specularReflectance">
			<xsl:choose>
				<xsl:when test="@type='microfacet'">1.0</xsl:when>
				<xsl:when test="@type='phong'">0.2</xsl:when>
				<xsl:when test="@type='ward'">0.2</xsl:when>
			</xsl:choose>
		</xsl:variable>
		<xsl:variable name="diffuseReflectance">
			<xsl:choose>
				<xsl:when test="@type='microfacet'">0.0</xsl:when>
				<xsl:when test="@type='phong'">0.5</xsl:when>
				<xsl:when test="@type='ward'">0.5</xsl:when>
			</xsl:choose>
		</xsl:variable>

		<bsdf>
			<xsl:apply-templates select="@*"/>

			<xsl:choose>
				<xsl:when test="node()[@name='diffuseReflectance']">
					<xsl:apply-templates select="node()[@name='diffuseReflectance']" mode="scaled">
						<xsl:with-param name="scale" select="$diffuseAmount"/>
					</xsl:apply-templates>
				</xsl:when>
				<xsl:otherwise>
					<spectrum name="diffuseReflectance" value="{number($diffuseAmount)*number($diffuseReflectance)}"/>
				</xsl:otherwise>
			</xsl:choose>

			<xsl:choose>
				<xsl:when test="node()[@name='specularReflectance']">
					<xsl:apply-templates select="node()[@name='specularReflectance']" mode="scaled">
						<xsl:with-param name="scale" select="$specularAmount"/>
					</xsl:apply-templates>
				</xsl:when>
				<xsl:otherwise>
					<spectrum name="specularReflectance" value="{number($specularAmount)*number($specularReflectance)}"/>
				</xsl:otherwise>
			</xsl:choose>
			<xsl:apply-templates select="node()[@name!='specularReflectance' and @name!='diffuseReflectance' and @name!='diffuseAmount' and @name!='specularAmount']"/>
		</bsdf>
	</xsl:template>

	<xsl:template match="ref|texture|rgb|srgb|spectrum|blackbody" mode="scaled">
		<xsl:param name="scale"/>
		<xsl:choose>
			<xsl:when test="number($scale)!=1">
				<texture type="scale">
					<xsl:if test="@name">
						<xsl:attribute name="name">
							<xsl:value-of select="@name"/>
						</xsl:attribute>
					</xsl:if>
					<float name="scale" value="{$scale}"/>
					<xsl:choose>
						<xsl:when test="texture|ref">
							<xsl:copy>
								<xsl:attribute name="name">value</xsl:attribute>
								<xsl:apply-templates select="@*[local-name() != 'name']|node()"/>
							</xsl:copy>
						</xsl:when>
						<xsl:otherwise>
							<xsl:copy>
								<xsl:attribute name="name">value</xsl:attribute>
								<xsl:apply-templates select="@*[local-name() != 'name']|node()"/>
							</xsl:copy>
						</xsl:otherwise>
					</xsl:choose>
				</texture>
			</xsl:when>
			<xsl:otherwise>
				<xsl:copy>
					<xsl:apply-templates select="@*|node()"/>
				</xsl:copy>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:template>

	<!-- Update the name of the lambertian plugin -->
	<xsl:template match="bsdf[@type='lambertian']/@type">
		<xsl:attribute name="type">diffuse</xsl:attribute>
	</xsl:template>
	
	<!-- Update the name of the microfacet plugin -->
	<xsl:template match="bsdf[@type='microfacet']/@type">
		<xsl:attribute name="type">roughplastic</xsl:attribute>
	</xsl:template>
	<xsl:template match="bsdf[@type='microfacet']/float[@name='alphaB']/@name">
		<xsl:attribute name="name">alpha</xsl:attribute>
	</xsl:template>
	
	<!-- There is no more 'mirror' plugin; replace with smooth chrome -->
	<xsl:template match="bsdf[@type='mirror']/@type">
		<xsl:attribute name="type">conductor</xsl:attribute>
	</xsl:template>
	<xsl:template match="bsdf[@type='mirror']">
		<xsl:copy>
			<xsl:apply-templates select="@*|node()"/>
			<string name="material" value="Cr"/>
		</xsl:copy>
	</xsl:template>

	<!-- Update the name of the roughmetal plugin -->
	<xsl:template match="bsdf[@type='roughmetal']/@type">
		<xsl:attribute name="type">roughconductor</xsl:attribute>
	</xsl:template>
	<xsl:template match="bsdf[@type='roughmetal']/float[@name='alphaB']/@name">
		<xsl:attribute name="name">alpha</xsl:attribute>
	</xsl:template>
	<xsl:template match="bsdf[@type='roughmetal']/float[@name='ior']/@name">
		<xsl:attribute name="name">eta</xsl:attribute>
	</xsl:template>

	<!-- Update the name of the roughglass plugin -->
	<xsl:template match="bsdf[@type='roughglass']/@type">
		<xsl:attribute name="type">roughdielectric</xsl:attribute>
	</xsl:template>
	<xsl:template match="bsdf[@type='roughglass']/float[@name='alphaB']/@name">
		<xsl:attribute name="name">alpha</xsl:attribute>
	</xsl:template>

	<!-- Update the name of the composite plugin -->
	<xsl:template match="bsdf[@type='composite']/@type">
		<xsl:attribute name="type">mixture</xsl:attribute>
	</xsl:template>
	
	<!-- Update the name of the exrtexture plugin -->
	<xsl:template match="texture[@type='exrtexture']/@type">
		<xsl:attribute name="type">bitmap</xsl:attribute>
	</xsl:template>

	<!-- Update the name of the ldrtexture plugin -->
	<xsl:template match="texture[@type='ldrtexture']/@type">
		<xsl:attribute name="type">bitmap</xsl:attribute>
	</xsl:template>

	<!-- Default copy rule -->
	<xsl:template match="@*|node()">
		<xsl:copy>
			<xsl:apply-templates select="@*|node()"/>
		</xsl:copy>
	</xsl:template>
</xsl:stylesheet>
